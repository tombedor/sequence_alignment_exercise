from Bio import SeqIO
class SequenceAssembler(object):
    def __init__(self, input_file):
        self.raw_sequences = SeqIO.parse(open(input_file),'fasta')
        self.build_sequence_linked_list_and_set_root()

    def build_sequence_linked_list_and_set_root(self):
        self.sequence_linked_list = {}
        for sequence in self.raw_sequences:
            seq_string = str(sequence.seq)
            seq_length = len(seq_string)
            self.sequence_linked_list[sequence.id] = {
                    'string':seq_string, 
                    'next_seq_name': None,
                    'next_seq_overlap_idx': None
                    }

        unmatched_end_seq_names= self.sequence_linked_list.keys()
        for start_name, start_seq in self.sequence_linked_list.iteritems():
            for end_name, end_seq in self.sequence_linked_list.iteritems():
                if start_name != end_name and end_name in unmatched_end_seq_names:
                    overlap_idx = self.get_overlap_idx(start_seq, end_seq)
                    if overlap_idx != -1:
                        self.sequence_linked_list[start_name]['next_seq_name'] = end_name
                        self.sequence_linked_list[start_name]['next_seq_overlap_idx'] = overlap_idx
                        unmatched_end_seq_names.remove(end_name)
                        break
        
        assert(len(unmatched_end_seq_names) == 1)
        self.root_name =unmatched_end_seq_names[0]
        return True

    def get_overlap_idx(self, start_seq, end_seq):
        start_seq_string = start_seq['string']
        end_seq_string = end_seq['string']
        min_overlap_length = max(len(start_seq_string), len(end_seq_string)) / 2
        min_overlap_string = end_seq_string[:min_overlap_length]

        overlap_index = start_seq_string.find(min_overlap_string)
        if overlap_index != -1 and len(start_seq_string[overlap_index:]) != len(min_overlap_string):
            full_overlap_string = end_seq_string[:overlap_index]
            return start_seq_string.find(full_overlap_string)
        else:
            return overlap_index


    def super_sequence(self):
        current_name = self.root_name
        super_seq = ""

        while current_name:
            current_seq = self.sequence_linked_list[current_name]
            truncate_idx_for_current_string = current_seq['next_seq_overlap_idx']
            if truncate_idx_for_current_string: 
                super_seq += current_seq['string'][:truncate_idx_for_current_string]
            else:
                super_seq += current_seq['string']
            current_name = current_seq.get('next_seq_name')
        return super_seq


class SuffixArray(object):
    def __init__(self, input_string):
        self.input_string = input_string
        self.array = range(len(input_string))
        self.array.sort(key=lambda i: input_string[i:])

    def suffix(self, array_index):
        string_idx = self.array[array_index]
        return self.input_string[string_idx:]

    def index(self, query_string, start_bound_array_idx = 0, end_bound_array_idx = None):
        if end_bound_array_idx is None:
            end_bound_array_idx = len(self.array) - 1

        search_idx = (end_bound_array_idx + start_bound_array_idx) / 2
        comparison_string = self.suffix(search_idx)

        if comparison_string == query_string:
            return self.array[search_idx]
        elif start_bound_array_idx == end_bound_array_idx:
            return -1
        elif comparison_string < query_string:
            return self.index(query_string, search_idx + 1, end_bound_array_idx)
        elif comparison_string > query_string:
            return self.index(query_string, start_bound_array_idx, search_idx - 1)
