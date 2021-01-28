# A Summary of the Mathematica Files

Running 11.2.0.0

* Type 1
  - groupmodk.nb: This file is supplementary material for proving what the commutator subgroup of Gamma(2)
  - tracesmodk.nb: This file is supplementary material for the proof for the set of traces modulo k in the commutator subgroup of Gamma(2)
  - verifiedCalculations: This file is supplementary material for chapter's 5 and 6 where calculations in the thesis are alluded to
  
* Type 2
  - class_number_search.nb: this search allows you to look for all of the conjugacy classes for given traces with smallest length. You are returned a csv that has all of the conjugacy class reps, their length, their decomposition, and their final position.
  - class_number_commutator.nb: this search allows you to look for all of the conjugacy classes for given traces with smallest length that lie in H'. You are returned a csv that has all of the conjugacy class reps, their length, their decomposition, and their final position.
  - existence_search.nb: this searches only the admissible traces to see if there exists at least one element in H' with that trace and if there doesn't it prints the trace that is a failure. So far we know 4.
  - where_python_left_off.nb: Uses properties of Markoff equations to narrow down the list given from python even more in deciding if something is a 1 or a 2 commutator.
  - homologyandlength.nb: Takes the representatives found in class_number_search.nb and exports a list of the homology class and length that the h conjugacy classes lie in.



