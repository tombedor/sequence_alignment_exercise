from Bio import SeqIO

input_file = './coding_challenge_data_set.txt'

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
sequence = fasta_sequences.next()
name = sequence.id
string = str(sequence.seq)
#SeqIO.to_dict?




