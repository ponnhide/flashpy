cimport cython
from cpython cimport array
from libcpp.vector cimport vector

@cython.boundscheck(False)
@cython.wraparound(False)
cdef hamming(array.array[char] x, array.array[char] y, int start1, int overlap, int start2):
    cdef int m 
    cdef double identity = 0.0
    for m in range(overlap): 
        if x[start1 + m] == y[start2 + m]:
            identity += 1.0
        else:
            pass
    identity = identity / overlap     
    return identity

@cython.boundscheck(False)
@cython.wraparound(False)
def merge(char *seq1, char *seq2, vector[int] score1, vector[int] score2, int min_overlap=50, int max_overlap=300, int allow_outies=1, double min_identity=0.5, double max_identity=1.0):
    ##Parameters
    cdef int reverse = 0 
    cdef int slide   = 0 
    cdef int overlap_length = 0
    cdef int start1 
    cdef int start2
    cdef array.array[char] sseq1, sseq2, template = array.array('b')
    cdef double identity
    cdef int i, j
    cdef int ls1 = len(seq1) 
    cdef int ls2 = len(seq2)
    cdef int current_slide   = 0
    cdef int current_overlap = 0
    cdef int current_start1  = 0
    cdef int current_start2  = 0
    cdef int left_length
    cdef int right_length
    cdef int outie
    cdef double current_identity, score_avg, cscore_avg
    cdef int  n = 0
    cdef vector[int] merged_score
    cdef array.array[char] merged_seq

    sseq1 = array.clone(template, ls1, False)
    for i in range(ls1):
        sseq1[i] = seq1[i]

    sseq2 = array.clone(template, ls2, False)
    for i in range(ls2): 
        sseq2[i] = seq2[i]

    current_overlap  = 0
    current_score    = None
    current_identity = min_identity
    current_slide    = 0 
    for i in range(ls1+ls2): 
        slide = i 
        if i < ls2 and allow_outies == 1:
            case = 0
            overlap_length = i 
            start1 = 0 
            start2 = ls2 - overlap_length
        elif i < ls1:
            case = 1
            overlap_length = ls2
            start1 = i 
            start2 = 0
        else:
            case = 2
            overlap_length = ls1 + ls2 - i 
            start1 = ls1 - overlap_length
            start2 = 0
        
        if min_overlap <= overlap_length <= max_overlap:
            identity = hamming(sseq1, sseq2, start1, overlap_length, start2)              
            if identity > current_identity or (identity == current_identity and current_overlap == 0):  
                current_slide     = slide   
                current_identity  = identity
                current_overlap   = overlap_length
                current_start1    = start1 
                current_start2    = start2

            elif identity == current_identity:
                n = 0
                cscore_avg = 0 
                for j in range(current_overlap):
                    if sseq1[current_start1+j] != sseq2[current_start2+j]:  
                        cscore_avg += max(score1[current_start1+j], score2[current_start2+j]) 
                        n += 1 
                csocre_avg = cscore_avg / n 
                
                n = 0
                score_avg = 0
                for j in range(overlap_length):
                    if sseq1[start1+j] != sseq2[start2+j]:
                        score_avg += max(score1[start1+j], score2[start2+j])
                        n += 1
                score_avg = score_avg / n 
                
                if score_avg > cscore_avg:
                    current_slide     = slide
                    current_identity  = identity
                    current_overlap   = overlap_length
                    current_start1    = start1 
                    current_start2    = start2
            else:
                pass
            
            if identity >= max_identity:
                break
    
    if current_identity < min_identity:
        return False
    else: 
        overlap_seq   = array.clone(array.array('b'), current_overlap, False)
        overlap_score = list(range(0, current_overlap)) 
        for i in range(current_overlap):  
            cn1 = score1[current_start1+i] 
            cn2 = score2[current_start2+i]
            cs1 = sseq1[current_start1+i]
            cs2 = sseq2[current_start2+i]
            if cs1 > cs2:
                overlap_score[i] = cn1
                overlap_seq[i]   = cs1 
            else:
                overlap_score[i] = cn2
                overlap_seq[i]   = cs2

        if current_slide < ls2 and allow_outies == True:
            outie = 1
            left_length  = ls2 - current_overlap
            right_length = ls1 - current_overlap 
            merged_seq   = array.clone(template, left_length + current_overlap + right_length, False) 
            merged_score = list(range(0, left_length + current_overlap + right_length)) 
            for i in range(left_length):
                merged_seq[i]   = sseq2[i]
                merged_score[i] = score2[i] 
            
            for i in range(current_overlap):
                merged_seq[i+left_length]   = overlap_seq[i]
                merged_score[i+left_length] = overlap_score[i]
            
            for i in range(right_length):
                merged_seq[i+left_length+current_overlap]   = sseq1[current_overlap + i]    
                merged_score[i+left_length+current_overlap] = score1[current_overlap + i]    
       
        else:
            outie = 0 
            left_length  = ls1 - current_overlap
            right_length = ls2 - current_overlap
            merged_seq   = array.clone(template, left_length + current_overlap + right_length, False) 
            merged_score = list(range(0, left_length + current_overlap + right_length))
            for i in range(left_length):
                merged_seq[i]   = sseq1[i]
                merged_score[i] = score1[i] 
            
            for i in range(current_overlap):
                merged_seq[i+left_length]   = overlap_seq[i]
                merged_score[i+left_length] = overlap_score[i]
            
            for i in range(right_length):
                merged_seq[i+left_length+current_overlap]   = sseq2[current_overlap + i]    
                merged_score[i+left_length+current_overlap] = score2[current_overlap + i]    
       
    return merged_seq, merged_score, outie, current_overlap, current_identity
