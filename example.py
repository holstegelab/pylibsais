import sys
import pylibsais
import numpy as np
from Bio import SeqIO


def parse_fasta(filename):
    """Parse a FASTA file and return a dictionary of sequences"""
    sequences = {}
    for record in SeqIO.parse(filename, 'fasta'):
        if record.seq != "":
            sequences[record.description] = str(record.seq)
    return sequences




def prepare_suffix_string(sequence_dict):
    """Prepare a sequence dictionary in format required by pylibsais."""
    keys = []
    values = []
    index = []
    pos = 0
    for k,v in sequence_dict.items():
        keys.append(k)
        values.append(v)
        pos += len(v) + 1
        index.append(pos)
    index = np.array(index)
    seq = '$'.join(values) + '$'

    return (seq, index) 


def get_positions(suffix_ar, kmer_idx, kmer_cnt):
    """Get all (sorted) positions of a kmer in a sequence"""
    positions = suffix_ar[kmer_idx:kmer_idx+kmer_cnt]
    return np.sort(positions)

def get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len):
    """Get the sequence of a kmer"""
    #get location of one of the kmer copies in original sequencedd
    xloc = suffix_ar[kmer_idx]
    #get the kmer sequence
    return seq[xloc:xloc+kmer_len]

def select_best_kmers(k_min, k_max, seq, index, min_count=2, min_consecutive=2, min_consecutive_bp=6):
    """Select k-mers based on the amount of sequence masked.

    :param seq: the sequence to search, as prepared by prepare_suffix_string
    :param index: the index of the end of each sequence in seq, as prepared by prepare_suffix_string
    """

    #create suffix and LCP array
    suffix_ar, lcp_ar = pylibsais.sais(seq)
    #determine maximum length of valid suffixes at each position (should stop at $ and # symbol)
    mkmer_ar = pylibsais.max_suffix(seq)
    

    #get all kmers with min_count or more copies
    #returns only maximal kmers (e.g. if each AG can be expanded to AGG, only AGG is returned)
    kmers = list(pylibsais.kmer_count(seq, suffix_ar, lcp_ar, mkmer_ar, k_min, k_max, min_count))
    kmers.sort(key=lambda x: (x[0] * x[2]), reverse=True) #sort on length * count, so that kmers that mask longer sequences are first

    #print(f'KMER CANDIDATES: {len(kmers)}')
    #for kmer_len, kmer_idx, kmer_cnt in kmers:
    #    print(f"- {get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len)}: {kmer_cnt} copies of length {kmer_len}")
    
    res = []
    max_continuous_masked_bp = 0

    #walk across possible kmers
    for kmer_len, kmer_idx, kmer_cnt in kmers:
        #stop if we cannot improve on the current best
        if (kmer_cnt * kmer_len) < max_continuous_masked_bp:
            break


        #determine how much of the sequence is masked by this kmer
        total_masked, max_indiv_seq_count, max_consecutive_count = pylibsais.kmer_mask_potential(suffix_ar, mkmer_ar, index, kmer_len, kmer_idx, kmer_cnt)
        
        #do not report kmer if it is worse than the current best
        #only report kmers that occur at least 2 times in a single sequence
        if max_consecutive_count * kmer_len < max_continuous_masked_bp or \
            max_indiv_seq_count == 1 or \
            max_consecutive_count < min_consecutive or \
            max_consecutive_count * kmer_len < min_consecutive_bp:
            continue

        max_continuous_masked_bp = max_consecutive_count * kmer_len

        #get the kmer sequence
        kmer_s = get_kmer_sequence(seq, suffix_ar, kmer_idx, kmer_len)
        min_kmer = pylibsais.min_string(kmer_s)
        
        #get all positions of the kmer in 'seq' (can be overlapping)
        positions = get_positions(suffix_ar, kmer_idx, kmer_cnt)

        res.append({'kmer':kmer_s, 'min_kmer': min_kmer, 'suffix_cnt': kmer_cnt, 'total_masked': total_masked, 
                        'max_indiv_seq_count':max_indiv_seq_count, 'max_consecutive_masked':max_consecutive_count * kmer_len, 'pos':positions, 'idx':kmer_idx})
   
    #sort kmers on priority: max continuous masked, then max count in individual sequence, then total masked, then length, then alphabetically
    res.sort(key=lambda x: (x['max_consecutive_masked'], x['max_indiv_seq_count'], x['total_masked'], len(x['kmer']), x['kmer']), reverse=True)

    return (res, suffix_ar, mkmer_ar)



sequence_dict = parse_fasta(sys.argv[1])

seq, index= prepare_suffix_string(sequence_dict)

#kmers that are selected
selected_kmers = []

#positions that are masked (list of tuples of position and kmer)
marked_positions = []

#repeat until no (consequtive) kmers are found
while True:    
    #min_consecutive=2: select kmers that occur at least in stretches of 2 consequtive k-mers.
    #min_count=2 indicates that only kmers are reported that occur at least 
    #2 times in a single sequence (superfluous as min_consecutive=2)
    res, sa, mask = select_best_kmers(2, 10, seq, index, min_count=2, min_consecutive=2)
    if len(res) == 0:
        break
    
    selected = res[0]
    selected_kmers.append(selected) 

    print(f"SELECT KMER: {selected['kmer']}")
    for k,v in selected.items():
        if k != 'kmer':
            print(f"- {k}: {v}")
    
    print('MASKED:')
    print(pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '.')[0])
    print('\n' * 2)
    
    #mask sequence with # symbol. The '2' indicates that only stretches of at least 2 consecutive kmers are masked.
    seq, marked_pos = pylibsais.kmer_mask(seq, sa, mask, len(selected['kmer']), selected['idx'], selected['suffix_cnt'], 2, '#')
    marked_positions.extend([(e, selected['kmer']) for e in marked_pos])
    #seq = seq.replace(selected['kmer'], '#' * len(selected['kmer']))


#mask now also single copies of the found kmers
for selected in selected_kmers:
    print(f"MASK KMER: {selected['kmer']}")
    print('MASKED:')
    print(pylibsais.kmer_mask_simple(seq, selected['kmer'], '.'))
    print('\n' * 2)
    #mask sequence with # symbol
    seq, marked_pos = pylibsais.kmer_mask_simple(seq, selected['kmer'], '#')
    marked_positions.extend([(e, selected['kmer']) for e in marked_pos])


for s in selected_kmers:
    print(s['kmer'])

marked_positions.sort(key=lambda x:x[0])
#print(marked_positions)
