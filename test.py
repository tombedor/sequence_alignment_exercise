from sequence_assembler import SequenceAssembler
from sequence_assembler import SuffixArray

assert(SequenceAssembler('./test_data_set.txt').super_sequence() == 'ATTAGACCTGCCGGAATAC')
print("Tests pass")
