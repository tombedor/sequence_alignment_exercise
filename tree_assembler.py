from Bio import SeqIO
class SequenceMap(object):
    def __init__(self, input_file):
        self.raw_sequences = SeqIO.parse(open(input_file),'fasta')
        self.build_sequence_graph_and_set_root()

    def build_sequence_graph_and_set_root(self):
        sequence_graph = {}
        for sequence in self.raw_sequences:
            seq_string = str(sequence.seq)
            seq_length = len(seq_string)
            sequence_graph[sequence.id] = {
                    'string':seq_string, 
                    }
        root_candidates = sequence_graph.keys()
        for start_name, start_seq in sequence_graph.iteritems():
            for end_name, end_seq in sequence_graph.iteritems():
                if start_name != end_name:
                    overlap_idx = self.get_overlap_idx(start_seq, end_seq)
                    if overlap_idx != -1:
                        sequence_graph[start_name]['next_seq_name'] = end_name
                        sequence_graph[start_name]['next_seq_overlap_idx'] = overlap_idx
                        root_candidates.remove(end_name)
        self.sequence_graph = sequence_graph
        
        assert(len(root_candidates) == 1)
        self.root_name = root_candidates[0]
        return True

    def get_overlap_idx(self, start_seq, end_seq):
        start_seq_string = start_seq['string']
        end_seq_string = end_seq['string']
        overlap_length = max(len(start_seq_string), len(end_seq_string)) / 2

        overlap_string = end_seq_string[:overlap_length]
        return start_seq_string.find(overlap_string)

    def super_sequence(self):
        current_name = self.root_name
        super_seq = ""

        while current_name:
            current_seq = self.sequence_graph[current_name]
            truncate_idx = current_seq.get('next_seq_overlap_idx', len(current_seq['string']))
            super_seq += current_seq['string'][:truncate_idx]
            current_name = current_seq.get('next_seq_name')
        return super_seq
