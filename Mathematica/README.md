# A Summary of the Mathematica Files

Running 11.2.0.0

* Type 1
  - groupmodk.nb: This file is supplmentary material for proving what the commutator subgroup of Gamma(2)
  - tracesmodk.nb: This file is supplmentary material for the proof for the set of traces modulo k in the commutator subgroup of Gamma(2)
  
* Type 2
  - class_number_search.nb: this search allows you to look for all of the conjugacy classes for given traces with the choice of choosing the smallest walk or just the first one found and whether to look only for admissible traces or all possible traces congruent to 2 modulo 16.  This can go very fast or very slow depending on which you choose. You are returned a csv that has all of the conjugacy class reps, their length, their decomposition, and their final position.
  - commutator_search_only.nb: this searches only the admissible traces to see if there exists at least one and if there doesn't it prints the trace that is a failure. So far we know 4.
  - where_python_left_off.nb: Uses properties of Markoff equations to see narrow down the list given from python even more in deciding if something is a 1 or a 2 commutator.



