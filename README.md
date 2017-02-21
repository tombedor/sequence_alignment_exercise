# Coding Challenge

## Summary

The prompt stated that the collection of subsequences could be combined to form exactly one supersequence. For this 
exercise, I assumed that the latter portion of each subsequence overlapped the front of exactly one other subsequence, 
aside from the end fragment. This assumption works if the subsequences are random and long enough. With real genomic
data, the presence of long repeated sequences might break this assumption.

## Algorithm

The data structure I used was a linked list. For each pair, I tried to overlay a "start" sequence on an "end" sequence.
If I found a match, I recorded the index of overlap, and removed the "end" sequence from a list of unmatched
subsequences. The last sequence left from this list was the root of the linked list.

Then, starting from the root, I concatenated each string (truncated at the point of overlap) to form the supersequence. 

To test whether one "start" sequence overlapped with an "end" sequence, I took the first half of the "end" sequence and 
tested whether it was a substring of the "start" sequence. If it overlapped at character N, I compared the last N 
characters of the "start" sequence with the first N characters of the "end" sequence.
 
I estimated the run time of this approach to be `O(n*n*l)`, where "n" is the number of sequences, and "l" is the
approximate length of each sequences. The reasoning is that each string comparison is about `O(l)` in the worst case, when 
only the last character differs between the two strings.

## Alternative Approach

I looked into other data structures that could improve on this, and implemented a suffix tree for string comparison.
For this, I first stored the sort order of each possible suffix of the "start" sequence. Using this array, a search for
the "end" fragment would take O(log l) time with binary search. Since the "end" fragment wouldn't usually be an exact 
suffix, when doing a binary search I truncated the search string to match the length of the "end" fragment.

This approach ended up being much slower than using Python's native string comparison operator, even without the cost of
building the suffix tree. I'm not sure why.

## Usage

The script uses BioPython to parse FASTA. To install:

`easy_install -f http://biopython.org/DIST/ biopython`

The supersequence is returned by:

`SequenceAssembler('./FASTA_FILE_NAME.txt').super_sequence()`

For testing, I used the data from the prompt email, and used the script `generate_test_data.py` to generate a few more
test scripts.
 
## Files
 
My solution is in `sequence_assembler.py`. The file can run the larger data set with:
`python run.py`

Tests can be run with:
`python test.py`

The file `alternative_sequence_assembler_with_suffix_array.py` has the alternative implementation with suffix arrays.
