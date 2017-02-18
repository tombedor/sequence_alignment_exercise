from tree_assembler import SequenceMap
mapper = SequenceMap('./test_data_set.txt')
assert(mapper.super_sequence() == 'ATTAGACCTGCCGGAATAC')
print "Test passes"

