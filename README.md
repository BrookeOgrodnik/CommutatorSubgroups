# On The Local-Global Conjecture for Commutator Subgroups
This folder contains Mathematica and Python Code that is used in the research for my PhD Thesis on the Local-Global Conjuecture for Thin Groups

## Algorithms in this repository
* **class_number_searches**
  - File Type: Mathematica
  - Summary: Uses relationship between Binary Quadratic Forms and Hyperbolic matrices to find all conjugacy classes (w.r.t. SL2(Z)) of the matrices in the subgroup, L, of SL2(Z) (that is the set of matrices congruent to I or 5I modulo 8) that have trace t where t is congruent to 2 modulo 16 and makes sure to return the ones whose decomposition starts with a power of {{1,2},{0,1}} and chooses the shortest length.
  
* **commutator_search_only**
  - File Type: Mathematica
  - Summary: Searches to see if a given admissible trace is the trace of an element in the commutator subgroup.
 
 
## Files from the Algorithms
* **info_on_traces**
  - File Type: CSV
  - Summary: The output from the class_number_search
