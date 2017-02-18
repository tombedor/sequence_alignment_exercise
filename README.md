Parse FASTA sequences using Biopyton. To install:
easy_install -f http://biopython.org/DIST/ biopython


bowtie
MUMmer

For each approach, I'll assume each of n subsequences has a length that is the same order of magnitude: l


construct a bidirectional graph, and track when a non-match or a match has been found.



Naive approach - nested loop using native Python string comparators:
  For each sequence i:
    For each sequence j:
      Using str.index(str), If the first half of string i in string j, match and remove strings from pool.

To do this, we must keep track of the minimum overlap distance we can use.
The native Python "in" operator runs in O(l) time:
https://wiki.python.org/moin/TimeComplexity

With each round, we'll have one less string to compare, so we'll have O(n!) comparisons to do. Total runtime is:
O(l*n^2)


suffix tree with graph

O(l*log l * n + log l * n^2)



    
