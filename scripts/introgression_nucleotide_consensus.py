import os
import argparse
from Bio import SeqIO
import logging
from typing import Dict, List, Tuple, Set, Optional
import numpy as np


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_group_file(group_file: str) -> Dict[str, List[str]]:
    groups = {}
    try:
        with open(group_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    logger.warning(f"Ignore lines with incorrect formatting: {line}")
                    continue
                seq_id = parts[0]
                group_name = parts[1]
                if group_name not in groups:
                    groups[group_name] = []
                groups[group_name].append(seq_id)
        
        # Filter out groups with a sequence count less than or equal to 3
        filtered_groups = {group: seqs for group, seqs in groups.items() if len(seqs) > 3}
        
        # Log the filtered groups
        for group, seqs in groups.items():
            if len(seqs) <= 3:
                logger.info(f"Ignored '{group}': sequence number {len(seqs)}，less than or equal to 3")
        
        logger.info(f"Read {len(groups)} groups from data file，retain {len(filtered_groups)} groups")
        for group_name, seq_ids in filtered_groups.items():
            logger.info(f"  group '{group_name}': {len(seq_ids)} sequences")
        
        return filtered_groups
    except Exception as e:
        logger.error(f"Error in reading data file: {e}")
        return {}

def extract_sequences(fasta_file: str) -> Dict[str, str]:
    sequences = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id
            if seq_id in sequences:
                logger.warning(f"Sequence ID '{seq_id}' is duplicated and will be overwritten")
            sequences[seq_id] = str(record.seq)
        logger.info(f"Read {len(sequences)} sequences from the FASTA file")
        return sequences
    except Exception as e:
        logger.error(f"Error in reading the FASTA file: {e}")
        return sequences

def generate_target_source_pairs(
    groups: Dict[str, List[str]]
) -> List[Tuple[str, str]]:
    pairs = []
    
    if len(groups) < 2:
        logger.error("There are not enough groups (at least 2 groups are required, and each group must have more than 3 sequences)")
        return []
    
    # For each group, treat it as a target group
    for tar_group in groups.keys():
        # Each of the remaining groups serves as a source group in turn
        for sour_group in groups.keys():
            if sour_group != tar_group:
                pairs.append((tar_group, sour_group))
    
    logger.info(f"generating {len(pairs)} source-target pairs")
    return pairs

def sliding_windows(sequence: str, window_size: int = 100) -> List[Tuple[int, int, str]]:
    windows = []
    seq_len = len(sequence)
    for start in range(0, seq_len, window_size):
        end = min(start + window_size, seq_len)
        window = sequence[start:end]
        if len(window) == window_size:
            windows.append((start, end, window))
    return windows

def calculate_nucleotide_identity(seq1: str, seq2: str) -> float:
    if len(seq1) != len(seq2):
        raise ValueError("The lengths of the sequences being compared are inconsistent")
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / len(seq1) * 100

def process_genome_pair(
    target_genome: Dict[str, str],
    source_genome: Dict[str, str],
    window_size: int = 100
) -> Dict[str, Dict[int, Dict[str, List[float]]]]:

    results = {}
    
    # For each sequence in target group
    for tar_seq_id, tar_seq in target_genome.items():
        results[tar_seq_id] = {}
        
        # Generating non-overlapping sliding windows
        windows = sliding_windows(tar_seq, window_size)
        
        # Calculating the consistency with all sequences in the source group for each window
        for start, end, tar_window_seq in windows:
            window_pos = start // window_size  
            results[tar_seq_id][window_pos] = {}
            
            for sour_seq_id, sour_seq in source_genome.items():
                sour_window_seq = sour_seq[start:end]
                
                if len(sour_window_seq) == window_size:
                    identity = calculate_nucleotide_identity(tar_window_seq, sour_window_seq)
                    
                    if sour_seq_id not in results[tar_seq_id][window_pos]:
                        results[tar_seq_id][window_pos][sour_seq_id] = []
                    results[tar_seq_id][window_pos][sour_seq_id].append(identity)
    
    return results

def calculate_target_internal_identity(
    target_sequences: Dict[str, str],
    window_size: int = 100
) -> Dict[int, Dict[str, List[float]]]:
    internal_identity = {}
    
    seq_ids = list(target_sequences.keys())
    
    for i in range(len(seq_ids)):
        for j in range(i+1, len(seq_ids)):
            seq_id1 = seq_ids[i]
            seq_id2 = seq_ids[j]
            
            seq1 = target_sequences[seq_id1]
            seq2 = target_sequences[seq_id2]
            
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]
            
            windows = sliding_windows(seq1, window_size)
            
            for start, end, _ in windows:
                window_pos = start // window_size
                
                window_seq1 = seq1[start:end]
                window_seq2 = seq2[start:end]
                
                identity = calculate_nucleotide_identity(window_seq1, window_seq2)
                
                if window_pos not in internal_identity:
                    internal_identity[window_pos] = {}
                
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

    window_stats = {}
    
    for window_pos, seq_identities in window_identities.items():
        all_identities = []
        for identities in seq_identities.values():
            all_identities.extend(identities)
        
        if all_identities:
            min_identity = min(all_identities)
            max_identity = max(all_identities)
            avg_identity = np.mean(all_identities)
            window_stats[window_pos] = (min_identity, max_identity, avg_identity)
    
    return window_stats

def detect_introgression_fragments(
    tar_sour_identity: Dict[str, Dict[int, Dict[str, List[float]]]],
    internal_identity: Dict[int, Dict[str, List[float]]],
    threshold: int,
    window_size: int = 100
) -> Set[int]:

    introgression_fragments = set()
    
    all_window_positions = set()
    for tar_seq_id, windows in tar_sour_identity.items():
        all_window_positions.update(windows.keys())
    
    for window_pos in all_window_positions:
        if window_pos in internal_identity:
            found = False
            
            for tar_seq_id, windows in tar_sour_identity.items():
                if window_pos in windows:
                    if tar_seq_id in internal_identity[window_pos]:
                        tar_internal_identities = internal_identity[window_pos][tar_seq_id]
                        
                        tar_internal_max = max(tar_internal_identities)
                        
                        sour_seqs = windows[window_pos]
                        tar_sour_identities = [ident for seqs in sour_seqs.values() for ident in seqs]
                        if tar_sour_identities:
                            tar_sour_max = max(tar_sour_identities)
                            
                            if tar_sour_max > tar_internal_max and tar_sour_max >= threshold:
                                introgression_fragments.add(window_pos)
                                found = True
                                break  
            
            if found:
                continue
    
    return introgression_fragments

def calculate_introgression_fraction(
    introgression_fragments: Set[int],
    total_windows: int
) -> float:
    if total_windows == 0:
        return 0.0
    
    return len(introgression_fragments) / total_windows * 100

def analyze_pair(
    sequences: Dict[str, str],
    groups: Dict[str, List[str]],
    tar_group: str,
    sour_group: str,
    output_dir: str,
    threshold: int,
    window_size: int = 100
) -> None:
    try:
        pair_dir = os.path.join(output_dir, f"{tar_group}_vs_{sour_group}")
        os.makedirs(pair_dir, exist_ok=True)
        
        tar_sequences = {seq_id: sequences[seq_id] for seq_id in groups[tar_group] if seq_id in sequences}
        sour_sequences = {seq_id: sequences[seq_id] for seq_id in groups[sour_group] if seq_id in sequences}
        
        logger.info(f"analyzing: {tar_group} (target) vs {sour_group} (source)")
        logger.info(f"  target sequences: {len(tar_sequences)} ")
        logger.info(f"  source sequences: {len(sour_sequences)} ")
        
        tar_fasta = os.path.join(pair_dir, "target.fasta")
        sour_fasta = os.path.join(pair_dir, "source.fasta")
        
        with open(tar_fasta, 'w') as f:
            for seq_id, seq in tar_sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        with open(sour_fasta, 'w') as f:
            for seq_id, seq in sour_sequences.items():
                f.write(f">{seq_id}\n{seq}\n")
        
        internal_identity = calculate_target_internal_identity(tar_sequences, window_size)
        
        internal_stats = calculate_window_statistics(internal_identity)
        
        tar_sour_identity = process_genome_pair(tar_sequences, sour_sequences, window_size)
        
        introgression_fragments = detect_introgression_fragments(
            tar_sour_identity, internal_identity, window_size
        )
        
        total_windows = 0
        for tar_seq_id, windows in tar_sour_identity.items():
            total_windows = len(windows)
        
        introgression_fraction = calculate_introgression_fraction(introgression_fragments, total_windows)
        logger.info(f"introgression_fraction (S_i): {introgression_fraction:.2f}%")
        
        detailed_results_file = os.path.join(pair_dir, "detailed_nucleotide_identity.csv")
        with open(detailed_results_file, 'w') as f:
            f.write("tar_Seq_ID,Window_Start,Window_End,sour_Seq_ID,Identity\n")
            
            for tar_seq_id, windows in tar_sour_identity.items():
                for window_pos, sour_seqs in windows.items():
                    start = window_pos * window_size
                    end = start + window_size
                    
                    for sour_seq_id, identities in sour_seqs.items():
                        for identity in identities:
                            f.write(f"{tar_seq_id},{start},{end},{sour_seq_id},{identity:.2f}\n")
        
        stats_results_file = os.path.join(pair_dir, "window_statistics.csv")
        with open(stats_results_file, 'w') as f:
            f.write("Window_Start,Window_End,tar_sour_Min,tar_sour_Max,tar_sour_Avg,")
            f.write("tar_Internal_Min,tar_Internal_Max,tar_Internal_Avg,Is_Introgressed\n")
            
            all_window_positions = set()
            for tar_seq_id, windows in tar_sour_identity.items():
                all_window_positions.update(windows.keys())
            
            for window_pos in sorted(all_window_positions):
                start = window_pos * window_size
                end = start + window_size
                
                all_tar_sour_identities = []
                for tar_seq_id, windows in tar_sour_identity.items():
                    if window_pos in windows:
                        sour_seqs = windows[window_pos]
                        tar_sour_identities = [ident for seqs in sour_seqs.values() for ident in seqs]
                        all_tar_sour_identities.extend(tar_sour_identities)
                
                if all_tar_sour_identities:
                    tar_sour_min = min(all_tar_sour_identities)
                    tar_sour_max = max(all_tar_sour_identities)
                    tar_sour_avg = np.mean(all_tar_sour_identities)
                else:
                    tar_sour_min, tar_sour_max, tar_sour_avg = float('nan'), float('nan'), float('nan')
                
                if window_pos in internal_stats:
                    tar_internal_min, tar_internal_max, tar_internal_avg = internal_stats[window_pos]
                else:
                    tar_internal_min, tar_internal_max, tar_internal_avg = float('nan'), float('nan'), float('nan')
                
                is_introgressed = window_pos in introgression_fragments
                
                f.write(f"{start},{end},{tar_sour_min:.2f},{tar_sour_max:.2f},{tar_sour_avg:.2f},")
                f.write(f"{tar_internal_min:.2f},{tar_internal_max:.2f},{tar_internal_avg:.2f},")
                f.write(f"{is_introgressed}\n")
        
        summary_file = os.path.join(pair_dir, "summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Target group\t{tar_group}\n")
            f.write(f"Source group\t{sour_group}\n")
            f.write(f"Number of target sequences\t{len(tar_sequences)}\n")
            f.write(f"Number of source sequences\t{len(sour_sequences)}\n")
            f.write(f"Total windows\t{total_windows}\n")
            f.write(f"Number of introgressed windows\t{len(introgression_fragments)}\n")
            f.write(f"Introgression fraction (Si)\t{introgression_fraction:.2f}%\n")
        
        logger.info(f"Finished, results were saved in: {pair_dir}")
        
    except Exception as e:
        logger.error(f"Error in {tar_group} vs {sour_group} : {e}")

def main():
    parser = argparse.ArgumentParser(description='Dynamically define reference and candidate species based on the grouping file and perform analysis')
    parser.add_argument('--fasta', required=True, help='FASTA file containing multiple core genome sequences')
    parser.add_argument('--group', required=True, help='Grouping file defining the group of each sequence')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--window-size', type=int, default=100, help='Sliding window size (default is 100bp)')
    parser.add_argument('--threshold', type=int, required=True, help='Source similarity threshold (default is 100%)')
    args = parser.parse_args()
    
    # Create the output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Read the grouping file
    groups = read_group_file(args.group)
    if not groups:
        logger.error("No valid grouping information was read or the number of sequences in all groups is less than or equal to 3. Program exiting.")
        return
    
    # Extract all sequences
    sequences = extract_sequences(args.fasta)
    if not sequences:
        logger.error("No valid sequence information was read. Program exiting")
        return
    
    # Generate all target-source pairs
    pairs = generate_target_source_pairs(groups)
    if not pairs:
        logger.error("No valid target-source pairs were generated. Program exiting")
        return
    
    # Analyze each pair
    for tar_group, sour_group in pairs:
        analyze_pair(
            sequences,
            groups,
            tar_group,
            sour_group,
            args.output,
            args.window_size,
            args.threshold
        )
    
    logger.info(f"All analyses completed. Results are saved in: {args.output}")

if __name__ == "__main__":
    main()
