import os
import argparse
from Bio import SeqIO
import logging
from typing import Dict, List, Tuple, Set, Optional
import numpy as np

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_group_file(group_file: str) -> Dict[str, List[str]]:
    """读取分组文件，返回序列分组信息"""
    groups = {}
    try:
        with open(group_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    logger.warning(f"忽略格式不正确的行: {line}")
                    continue
                seq_id = parts[0]
                group_name = parts[1]
                if group_name not in groups:
                    groups[group_name] = []
                groups[group_name].append(seq_id)
        
        # 过滤掉序列数少于等于3的分组
        filtered_groups = {group: seqs for group, seqs in groups.items() if len(seqs) > 3}
        
        # 记录被过滤的分组
        for group, seqs in groups.items():
            if len(seqs) <= 3:
                logger.info(f"忽略分组 '{group}': 序列数为 {len(seqs)}，少于等于3")
        
        logger.info(f"从分组文件读取到 {len(groups)} 个分组，保留 {len(filtered_groups)} 个分组")
        for group_name, seq_ids in filtered_groups.items():
            logger.info(f"  分组 '{group_name}': {len(seq_ids)} 条序列")
        
        return filtered_groups
    except Exception as e:
        logger.error(f"读取分组文件时出错: {e}")
        return {}

def extract_sequences(fasta_file: str) -> Dict[str, str]:
    """从FASTA文件中提取所有序列"""
    sequences = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            if seq_id in sequences:
                logger.warning(f"序列ID '{seq_id}' 重复，将被覆盖")
            sequences[seq_id] = str(record.seq)
        logger.info(f"从FASTA文件中读取了 {len(sequences)} 条序列")
        return sequences
    except Exception as e:
        logger.error(f"读取FASTA文件时出错: {e}")
        return sequences

def generate_reference_candidate_pairs(
    groups: Dict[str, List[str]]
) -> List[Tuple[str, str]]:
    """生成所有可能的参考物种和候选物种组合"""
    pairs = []
    
    # 检查是否有足够的分组
    if len(groups) < 2:
        logger.error("没有足够的分组（至少需要2个分组，且每个分组序列数大于3）")
        return []
    
    # 对每个分组，将其作为参考物种
    for ref_group in groups.keys():
        # 剩下的每个分组依次作为候选物种
        for cand_group in groups.keys():
            if cand_group != ref_group:
                pairs.append((ref_group, cand_group))
    
    logger.info(f"生成了 {len(pairs)} 个参考-候选物种组合")
    return pairs

def sliding_windows(sequence: str, window_size: int = 100) -> List[Tuple[int, int, str]]:
    """生成非重叠的滑动窗口"""
    windows = []
    seq_len = len(sequence)
    for start in range(0, seq_len, window_size):
        end = min(start + window_size, seq_len)
        window = sequence[start:end]
        # 只包含完整的窗口（长度为window_size）
        if len(window) == window_size:
            windows.append((start, end, window))
    return windows

def calculate_nucleotide_identity(seq1: str, seq2: str) -> float:
    """计算两个序列的核苷酸一致性"""
    if len(seq1) != len(seq2):
        raise ValueError("比较的序列长度不一致")
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / len(seq1) * 100

def process_genome_pair(
    reference_genome: Dict[str, str],
    candidate_genome: Dict[str, str],
    window_size: int = 100
) -> Dict[str, Dict[int, Dict[str, List[float]]]]:
    """
    计算参考基因组中每个序列与候选基因组中所有序列的一致性
    返回每个窗口的详细一致性数据
    """
    results = {}
    
    # 对参考基因组的每个序列
    for ref_seq_id, ref_seq in reference_genome.items():
        results[ref_seq_id] = {}
        
        # 生成滑动窗口
        windows = sliding_windows(ref_seq, window_size)
        
        # 对每个窗口，计算与候选基因组中所有序列的一致性
        for start, end, ref_window_seq in windows:
            window_pos = start // window_size  # 窗口位置索引
            results[ref_seq_id][window_pos] = {}
            
            # 计算与候选基因组中所有序列的一致性
            for cand_seq_id, cand_seq in candidate_genome.items():
                # 提取候选序列的对应窗口
                cand_window_seq = cand_seq[start:end]
                
                # 只处理长度匹配的窗口
                if len(cand_window_seq) == window_size:
                    # 计算一致性
                    identity = calculate_nucleotide_identity(ref_window_seq, cand_window_seq)
                    
                    # 存储结果
                    if cand_seq_id not in results[ref_seq_id][window_pos]:
                        results[ref_seq_id][window_pos][cand_seq_id] = []
                    results[ref_seq_id][window_pos][cand_seq_id].append(identity)
    
    return results

def calculate_reference_internal_identity(
    reference_sequences: Dict[str, str],
    window_size: int = 100
) -> Dict[int, Dict[str, List[float]]]:
    """
    计算参考基因组内部的核苷酸一致性，返回每个窗口每个参考序列与其他参考序列的一致性列表
    """
    internal_identity = {}
    
    # 获取参考基因组中的所有序列ID
    seq_ids = list(reference_sequences.keys())
    
    # 如果参考基因组中序列数量不足2条，无法计算内部一致性
    if len(seq_ids) < 2:
        logger.warning("参考基因组中序列数量不足，无法计算内部一致性")
        return internal_identity
    
    # 计算所有序列对之间的一致性
    for i in range(len(seq_ids)):
        for j in range(i+1, len(seq_ids)):
            seq_id1 = seq_ids[i]
            seq_id2 = seq_ids[j]
            
            seq1 = reference_sequences[seq_id1]
            seq2 = reference_sequences[seq_id2]
            
            # 确保两条序列长度相同
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]
            
            # 生成滑动窗口
            windows = sliding_windows(seq1, window_size)
            
            for start, end, _ in windows:
                window_pos = start // window_size
                
                # 提取窗口序列
                window_seq1 = seq1[start:end]
                window_seq2 = seq2[start:end]
                
                # 计算一致性
                identity = calculate_nucleotide_identity(window_seq1, window_seq2)
                
                # 更新窗口一致性结果
                if window_pos not in internal_identity:
                    internal_identity[window_pos] = {}
                
                # 存储序列对的一致性
                if seq_id1 not in internal_identity[window_pos]:
                    internal_identity[window_pos][seq_id1] = []
                if seq_id2 not in internal_identity[window_pos]:
                    internal_identity[window_pos][seq_id2] = []
                
                internal_identity[window_pos][seq_id1].append(identity)
                internal_identity[window_pos][seq_id2].append(identity)
    
    return internal_identity

def calculate_window_statistics(
    window_identities: Dict[int, Dict[str, List[float]]]
) -> Dict[int, Tuple[float, float, float]]:
    """
    计算每个窗口的统计量（最小、最大、平均一致性）
    """
    window_stats = {}
    
    for window_pos, seq_identities in window_identities.items():
        # 收集该窗口的所有一致性值
        all_identities = []
        for identities in seq_identities.values():
            all_identities.extend(identities)
        
        # 计算统计量
        if all_identities:
            min_identity = min(all_identities)
            max_identity = max(all_identities)
            avg_identity = np.mean(all_identities)
            window_stats[window_pos] = (min_identity, max_identity, avg_identity)
    
    return window_stats

def detect_introgression_fragments(
    ref_cand_identity: Dict[str, Dict[int, Dict[str, List[float]]]],
    internal_identity: Dict[int, Dict[str, List[float]]],
    threshold: int,
    window_size: int = 100
) -> Set[int]:
    """
    检测渐渗片段（当至少一个参考序列与候选序列的一致性最大值高于该参考序列内部一致性的最大值时认为发生渐渗）
    优化版本：只需要找到一个满足条件的参考序列即可将窗口标记为渐渗，无需存储所有参考序列的结果
    """
    introgression_fragments = set()
    
    # 遍历所有窗口位置
    all_window_positions = set()
    for ref_seq_id, windows in ref_cand_identity.items():
        all_window_positions.update(windows.keys())
    
    for window_pos in all_window_positions:
        # 检查该窗口是否有参考内部一致性数据
        if window_pos in internal_identity:
            # 标记是否找到满足条件的参考序列
            found = False
            
            # 遍历所有参考序列
            for ref_seq_id, windows in ref_cand_identity.items():
                # 检查该参考序列在当前窗口是否有候选-参考一致性数据
                if window_pos in windows:
                    # 获取参考序列在该窗口的内部一致性列表
                    if ref_seq_id in internal_identity[window_pos]:
                        ref_internal_identities = internal_identity[window_pos][ref_seq_id]
                        
                        # 计算参考序列内部一致性的最大值
                        ref_internal_max = max(ref_internal_identities)
                        
                        # 计算参考序列与候选序列的一致性的最大值
                        cand_seqs = windows[window_pos]
                        ref_cand_identities = [ident for seqs in cand_seqs.values() for ident in seqs]
                        if ref_cand_identities:
                            ref_cand_max = max(ref_cand_identities)
                            
                            # 如果参考-候选一致性最大值大于参考内部一致性最大值，则标记该窗口为渐渗片段
                            if ref_cand_max > ref_internal_max and ref_cand_max >= threshold:
                                introgression_fragments.add(window_pos)
                                found = True
                                break  # 找到一个满足条件的参考序列即可退出循环
            
            # 如果已经找到满足条件的参考序列，无需继续检查其他参考序列
            if found:
                continue
    
    return introgression_fragments

def calculate_introgression_fraction(
    introgression_fragments: Set[int],
    total_windows: int
) -> float:
    """计算渐渗分数"""
    if total_windows == 0:
        return 0.0
    
    return len(introgression_fragments) / total_windows * 100

def analyze_pair(
    sequences: Dict[str, str],
    groups: Dict[str, List[str]],
    ref_group: str,
    cand_group: str,
    output_dir: str,
    threshold: int,
    window_size: int = 100
) -> None:
    """分析一个参考-候选物种组合"""
    try:
        # 创建输出目录
        pair_dir = os.path.join(output_dir, f"{ref_group}_vs_{cand_group}")
        os.makedirs(pair_dir, exist_ok=True)
        
        # 获取参考和候选物种的序列
        ref_sequences = {seq_id: sequences[seq_id] for seq_id in groups[ref_group] if seq_id in sequences}
        cand_sequences = {seq_id: sequences[seq_id] for seq_id in groups[cand_group] if seq_id in sequences}
        
        logger.info(f"分析组合: {ref_group} (参考) vs {cand_group} (候选)")
        logger.info(f"  参考序列: {len(ref_sequences)} 条")
        logger.info(f"  候选序列: {len(cand_sequences)} 条")
        
        # 输出参考和候选序列到文件
        ref_fasta = os.path.join(pair_dir, "reference.fasta")
        cand_fasta = os.path.join(pair_dir, "candidate.fasta")
        
        with open(ref_fasta, 'w') as f:
            for seq_id, seq in ref_sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        with open(cand_fasta, 'w') as f:
            for seq_id, seq in cand_sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        # 计算参考基因组内部一致性
        logger.info("开始计算参考基因组内部一致性...")
        internal_identity = calculate_reference_internal_identity(ref_sequences, window_size)
        
        # 计算参考基因组内部一致性统计量
        internal_stats = calculate_window_statistics(internal_identity)
        
        # 处理基因组对，计算参考基因组中每个序列与候选基因组中所有序列的一致性
        logger.info("开始计算参考-候选基因组一致性...")
        ref_cand_identity = process_genome_pair(ref_sequences, cand_sequences, window_size)
        
        # 检测渐渗片段
        logger.info("开始检测渐渗片段...")
        introgression_fragments = detect_introgression_fragments(
            ref_cand_identity, internal_identity, window_size
        )
        
        # 计算总窗口数
        total_windows = 0
        for ref_seq_id, windows in ref_cand_identity.items():
            total_windows = len(windows)
        
        # 计算渐渗分数
        introgression_fraction = calculate_introgression_fraction(introgression_fragments, total_windows)
        logger.info(f"渐渗分数 (S_i): {introgression_fraction:.2f}%")
        
        # 输出详细的一致性结果
        detailed_results_file = os.path.join(pair_dir, "detailed_nucleotide_identity.csv")
        with open(detailed_results_file, 'w') as f:
            # 写入表头
            f.write("Ref_Seq_ID,Window_Start,Window_End,Cand_Seq_ID,Identity\n")
            
            # 写入详细的一致性数据
            for ref_seq_id, windows in ref_cand_identity.items():
                for window_pos, cand_seqs in windows.items():
                    start = window_pos * window_size
                    end = start + window_size
                    
                    for cand_seq_id, identities in cand_seqs.items():
                        # 每个窗口每个候选序列可能有多个一致性值
                        for identity in identities:
                            f.write(f"{ref_seq_id},{start},{end},{cand_seq_id},{identity:.2f}\n")
        
        # 输出窗口统计结果
        stats_results_file = os.path.join(pair_dir, "window_statistics.csv")
        with open(stats_results_file, 'w') as f:
            # 写入表头
            f.write("Window_Start,Window_End,Ref_Cand_Min,Ref_Cand_Max,Ref_Cand_Avg,")
            f.write("Ref_Internal_Min,Ref_Internal_Max,Ref_Internal_Avg,Is_Introgressed\n")
            
            # 写入窗口统计数据
            all_window_positions = set()
            for ref_seq_id, windows in ref_cand_identity.items():
                all_window_positions.update(windows.keys())
            
            for window_pos in sorted(all_window_positions):
                start = window_pos * window_size
                end = start + window_size
                
                # 计算所有参考序列在该窗口的参考-候选一致性统计量
                all_ref_cand_identities = []
                for ref_seq_id, windows in ref_cand_identity.items():
                    if window_pos in windows:
                        cand_seqs = windows[window_pos]
                        ref_cand_identities = [ident for seqs in cand_seqs.values() for ident in seqs]
                        all_ref_cand_identities.extend(ref_cand_identities)
                
                if all_ref_cand_identities:
                    ref_cand_min = min(all_ref_cand_identities)
                    ref_cand_max = max(all_ref_cand_identities)
                    ref_cand_avg = np.mean(all_ref_cand_identities)
                else:
                    ref_cand_min, ref_cand_max, ref_cand_avg = float('nan'), float('nan'), float('nan')
                
                # 获取参考内部一致性统计数据
                if window_pos in internal_stats:
                    ref_internal_min, ref_internal_max, ref_internal_avg = internal_stats[window_pos]
                else:
                    ref_internal_min, ref_internal_max, ref_internal_avg = float('nan'), float('nan'), float('nan')
                
                # 检查是否是渐渗片段
                is_introgressed = window_pos in introgression_fragments
                
                f.write(f"{start},{end},{ref_cand_min:.2f},{ref_cand_max:.2f},{ref_cand_avg:.2f},")
                f.write(f"{ref_internal_min:.2f},{ref_internal_max:.2f},{ref_internal_avg:.2f},")
                f.write(f"{is_introgressed}\n")
        
        # 输出渐渗分数汇总
        summary_file = os.path.join(pair_dir, "summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"参考物种\t{ref_group}\n")
            f.write(f"候选物种\t{cand_group}\n")
            f.write(f"参考序列数\t{len(ref_sequences)}\n")
            f.write(f"候选序列数\t{len(cand_sequences)}\n")
            f.write(f"总窗口数\t{total_windows}\n")
            f.write(f"渐渗窗口数\t{len(introgression_fragments)}\n")
            f.write(f"渐渗分数(Si)\t{introgression_fraction:.2f}%\n")
        
        logger.info(f"分析完成，结果保存在: {pair_dir}")
        
    except Exception as e:
        logger.error(f"分析组合 {ref_group} vs {cand_group} 时出错: {e}")

def main():
    parser = argparse.ArgumentParser(description='根据分组文件动态定义参考物种和候选物种并分析')
    parser.add_argument('--fasta', required=True, help='包含多条核心基因组序列的FASTA文件')
    parser.add_argument('--group', required=True, help='分组文件，定义每条序列的分组')
    parser.add_argument('--output', required=True, help='输出目录')
    parser.add_argument('--window-size', type=int, default=100, help='滑动窗口大小 (默认为100bp)')
    parser.add_argument('--threshold', type=int, required=True, help='候选相似度阈值 (默认为100%)')
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output, exist_ok=True)
    
    # 读取分组文件
    groups = read_group_file(args.group)
    if not groups:
        logger.error("没有读取到有效的分组信息或所有分组序列数均小于等于3，程序退出")
        return
    
    # 提取所有序列
    sequences = extract_sequences(args.fasta)
    if not sequences:
        logger.error("没有读取到有效的序列信息，程序退出")
        return
    
    # 生成所有参考-候选物种组合
    pairs = generate_reference_candidate_pairs(groups)
    if not pairs:
        logger.error("没有生成有效的参考-候选物种组合，程序退出")
        return
    
    # 分析每个组合
    for ref_group, cand_group in pairs:
        analyze_pair(
            sequences,
            groups,
            ref_group,
            cand_group,
            args.output,
            args.window_size,
            args.threshold
        )
    
    logger.info(f"所有分析完成，结果保存在: {args.output}")

if __name__ == "__main__":
    main()
