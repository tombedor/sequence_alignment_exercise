from tree_assembler import SequenceMap
from tree_assembler import SuffixArray

mapper = SequenceMap('./test_data_set.txt')
assert(mapper.super_sequence() == 'ATTAGACCTGCCGGAATAC')


suffix_array = SuffixArray('zadc')

assert(suffix_array.array == [1, 3, 2, 0])
assert(suffix_array.index('adc') == 1)



assert(suffix_array.index('ffff') == -1)
assert(suffix_array.index('c') == 3)



print "Test passes"

