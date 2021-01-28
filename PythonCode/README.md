# A Summary of the Python Files

Running python-3.7.6

* Type 1
  - Class_number.ipynb: Just takes the file trace_only_converted.csv and draws plots of the class number of L, the class number of the commutator subgroup, ratios, average lenght of walks, etc.
  - Walks_of_traces.ipynb: Takes all_traces_converted.csv and plots walks and heat maps for the given conjugacy classes
  - Walks_of_trace_hclass.ipynb: Does the same as Walks_of_traces.ipynb but does it for a file in a different format that represents the H conjugacy classes
  
  
* Type 2
  - Convert_files.ipynb: Takes the results of the mathematica code, named info_on_traces.csv and gets it into a form that works better in python
  - Genus_of_traces.ipynb: Does a preliminary search through admissible_trace_only_converted.csv to find the genus of the given traces.  If it cannot be concluded it is indicated in the output that is then sent to 1 last mathematica file where_python_left_off.nb. There is a search first that bounded each trace's genus above that should always be ran so you know that the options are only 1 or 2.
