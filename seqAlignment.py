from global_alignment import compute_global_alignment
from local_alignment import compute_local_alignment

sequence1 = 'VIVALASVEGAS'
sequence2 = 'VIVADAVIS'
gap = 1
match = 1
mismatch = 1

global_alignment = compute_global_alignment(sequence1, sequence2, gap, match, mismatch)
print(global_alignment)

print('----------------------------------------------------------')
sequence2 = 'TGTTACGG'
sequence1 = 'GGTTGACTA'
gap = 2
match = 3
mismatch = 3

local_alignment = compute_local_alignment(sequence1, sequence2, gap, match, mismatch)
print(local_alignment)
