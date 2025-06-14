import pandas as pd
from joblib import Parallel, delayed
from joblib import Memory
import time

def process_chunk(chunk, ps_values):
    try:
        mask = chunk['Site1'].isin(ps_values) & chunk['Site2'].isin(ps_values)
        filtered_chunk = chunk[mask]
        return filtered_chunk
    except Exception as e:
        print(f"error: {e}")
        return pd.DataFrame()



try:
    file1 = pd.read_table('sig_SNP_strict.txt')
    ps_values = set(file1['ps'])
except FileNotFoundError:
    print("file1 not found, please check")

chunk_size = 100000000  
results = []

chunks = pd.read_table('../snp.split.filteredMAF.gap5_split.LD.gz', chunksize=chunk_size)
start_time=time.time()
for i in range(1,10000):
    try:
        chunk=next(chunks)
        results.append(process_chunk(chunk, ps_values))
        del chunk
    except StopIteration:
         print(f"iteration times: {i}")
         break

end_time=time.time()
elapsed_time = end_time - start_time
print(f"elapsed time: {elapsed_time} 秒")

result = pd.concat(results, ignore_index=True)

output_file = 'strictSNP_LD_v2.txt'
result.to_csv(output_file, sep='\t', index=False)
print(f"result has been saved to {output_file}")