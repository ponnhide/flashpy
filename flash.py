import os 
import sys
import gzip 
import collections
import numpy as np
from tqdm import tqdm
import _flash as fl

hamming = lambda x,y: sum(tuple(map(lambda a,b: 1 if a==b else 0, x, y))) / len(x)
convert_ascii = lambda x: [ord(asc) - 33 for asc in x]
def read_fastq(fastq_name):
    """Read fastq file
    """
    seq_dict = collections.defaultdict(dict)
    if fastq_name.split(".")[-1] == "gz":
        f = gzip.open(fastq_name.replace("'","").replace("\\",""), mode="rt", encoding='utf-8')
    else:
        f = open(fastq_name.replace("'","").replace("\\",""))  
    n = 0 
    for line in f:
        if line[0] == "@" and n % 4 == 0:
            key = line[1:].rstrip() 
            key = key.split(" ")[0] 
            key = key.replace(":","_")
            seq_dict[key]["key"] = line[1:].rstrip() 
        elif n%4 == 1:
            seq_dict[key]["seq"] = line.rstrip()
        elif n%4 == 2:
            seq_dict[key]["option"] = line.rstrip() 
        elif n%4 == 3:
            seq_dict[key]["quality"] = [ord(asc) - 33 for asc in line.rstrip()] 
        n += 1
    f.close() 
    return seq_dict 

def merge(seq1, seq2, score1, score2, min_overlap=50, max_overlap=300, allow_outies=True, min_identity=0.5, max_identity=1.0, cython=True):
    #Pure python code to merge paired sequence. The execution speed is too slow, please use cythonized version fl.merge 
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    reverse = 0 
    if len(seq1) >= len(seq2):
        seq2 = seq2.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        score2 = score2[::-1]
    else:
        reverse = 1
        seq1 = seq1.translate(str.maketrans("ATGCRYKMSWBDHV","TACGYRMKWSVHDB"))[::-1]
        seq2, seq1 = seq1, seq2
    
    if cython == False:
        current_overlap  = 0
        current_score    = None
        current_identity = min_identity
        current_slide     = 0 
        for i in range(len(seq1)+len(seq2)): 
            slide = i 
            if i < len(seq2) and allow_outies == True:
                overlap_length = i 
                subseq1   = seq1[:overlap_length] 
                subseq2   = seq2[-1*overlap_length:]
                subscore1 = score1[:overlap_length]
                subscore2 = score2[-1*overlap_length:]   
            elif i < len(seq1):
                overlap_length = len(seq2) 
                subseq1   = seq1[i-len(seq2):i] 
                subseq2   = seq2
                subscore1 = score1[i-len(seq2):i]
                subscore2 = score2 
            else:
                overlap_length = len(seq1) + len(seq2) - i 
                subseq1   = seq1[-1*overlap_length:] 
                subseq2   = seq2[:overlap_length:]   
                subscore1 = score1[-1*overlap_length:] 
                subscore2 = score2[:overlap_length]
            
            if min_overlap <= overlap_length <= max_overlap:
                identity = hamming(subseq1, subseq2)             
                if identity > current_identity or (identity == current_identity and current_overlap == 0):
                    current_slide     = slide   
                    current_identity  = identity
                    current_overlap   = overlap_length
                    current_subseq1   = subseq1 
                    current_subseq2   = subseq2
                    current_subscore1 = subscore1 
                    current_subscore2 = subscore2
                    #for n1, n2, s1, s2 in zip(subseq1, subseq2, subscore1, subscore2):

                elif identity == current_identity:
                    n = 0
                    cscore_avg = 0 
                    for cn1, cn2, cs1, cs2 in zip(current_subseq1, current_subseq2, current_subscore1, current_subscore2):
                        if cn1 != cn2:  
                            cscore_avg += max(cs1, cs2) 
                            n += 1
                    csocre_avg = cscore_avg / n
                    
                    n = 0 
                    score_avg = 0
                    for n1, n2, s1, s2 in zip(subseq1, subseq2, subscore1, subscore2):
                        if n1 != n2:
                            score_sum += max(s1, s2)
                            n += 1
                    score_avg = score_avg / n 
                    if score_avg > cscore_avg:
                        current_slide     = slide
                        current_identity  = identity
                        current_overlap   = overlap_length
                        current_subseq1   = subseq1 
                        current_subseq2   = subseq2
                        current_subscore1 = subscore1 
                        current_subscore2 = subscore2
                else:
                    pass
                
                if identity >= max_identity:
                    break
        
        if current_identity < min_identity:
            return False
        else: 
            overlap_seq   = "" 
            overlap_score = []
            for cn1, cn2, cs1, cs2 in zip(current_subseq1, current_subseq2, current_subscore1, current_subscore2):
                if cs1 > cs2:
                    overlap_seq += cn1
                    overlap_score.append(cs1)
                else:
                    overlap_seq += cn2
                    overlap_score.append(cs2)
            
            if current_slide < len(seq2) and allow_outies == True:
                left_seq    = seq2[:-1*overlap_length]
                right_seq   = seq1[overlap_length:] 
                left_score  = score2[:-1*overlap_length]
                right_score = score1[overlap_length:] 
                        
            else:
                left_seq    = seq1[:-1*overlap_length]
                right_seq   = seq2[overlap_length:] 
                left_score  = score1[:-1*overlap_length]
                right_score = score2[overlap_length:] 
            merged_seq   = left_seq + overlap_seq + right_seq
            merged_score = left_score + overlap_score + right_score
    else:
        merged_seq, merged_score, current_slide, current_overlap, current_identity = fl.merge(seq1.encode('utf-8'), seq2.encode('utf-8'), score1, score2, min_overlap, max_overlap, allow_outies, min_identity, max_identity)
    return merged_seq, merged_score, current_slide, current_overlap, current_identity

def flash(read1, read2, min_overlap=10, max_overlap=300, allow_outies=True, min_identity=0.5, max_identity=1.0, show_progress=True, key_check=True):
    r1_dict = read_fastq(read1)
    r2_dict = read_fastq(read2) 
    
    dist_dict = collections.defaultdict(lambda:[0, 0])  
    merged_dict = collections.defaultdict(dict)
    if show_progress == True:
        if key_check == True:
            keys = [key for key in r1_dict.keys() if key in r2_dict]
            keys = tqdm(keys)
        else:
            keys = tqdm(r1_dict.keys(), total=len(r1_dict))
    else:
        if key_check == True:
            keys = [key for key in r1_dict.keys() if key in r2_dict]
        else:
            keys = r1_dict.keys()
    
    for key in keys:
        seq1   = r1_dict[key]["seq"]
        seq2   = r2_dict[key]["seq"]
        score1 = r1_dict[key]["quality"]
        score2 = r2_dict[key]["quality"]
        result = merge(seq1, seq2, score1, score2, min_overlap, max_overlap, allow_outies, min_identity, max_identity) 
        if result != False:
            merged_dict[key]["seq"]      = "".join(map(chr, result[0])) 
            merged_dict[key]["quality"]  = result[1]
            merged_dict[key]["r1_key"]   = r1_dict[key]["key"]
            merged_dict[key]["r2_key"]   = r2_dict[key]["key"]
            merged_dict[key]["identity"] = result[4]
            if result[2] == 1: 
                dist_dict["outie", result[3]][0] += 1
                dist_dict["outie", result[3]][1] += result[4]
            else:
                dist_dict["innie", result[3]][0] += 1
                dist_dict["innie", result[3]][1] += result[4]    
    for key in dist_dict:
        dist_dict[key][1] = dist_dict[key][1] / dist_dict[key][0] 
    
    return merged_dict, dist_dict

if __name__ == "__main__":
    #path = "../../sample_data/Target-ACE_sequence_data/Rep1/AID/polyC/chr19_8211972/" #R1.fastq, R2.fastq 
    merged_dict, dist_dict = flash("R1.fastq.gz", "R2.fastq.gz") 
    for key in dist_dict:
        print(key, *dist_dict[key], sep="\t") 
    
    n = 0
    for key in merged_dict:
        print(key, merged_dict[key]["seq"], sep="\t")
        n += 1
        if n == 10:
            break
