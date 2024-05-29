# File descriptions:  

* **Publication data**: Folder containing raw data used in our preprint (and submitted publication) on this work.
  
* **LHS_Sampling_And_Parameter_Export.ipynb**: This is a Jupyter notebook which conducts parameter sampling and generates .csv files containing, in each row of each file, a unique set of parameters to run in a simulation.  

* **Param_explore.py**: This is a Python script that iterates through all of the parameter sets contained in a single .csv file, runs the corresponding simulations, and outputs the results (along with a flag of whether or not the simulation reached steady state).  

* **Parameter_exploration.slurm**: This UNIX file distributes the individual param_explore.py jobs to the computing cluster. 

* **Stitcher.py**: This Python file collects all of the individual results files together and turns them into a single .csv.  

# Instructions 

1. The **LHS_Sampling_And_Parameter_Export.ipynb** is used to generate a list of parameter sets to run through the C. merolae compartmental model, based off of user-defined parameter upper and lower bounds and a Latin Hypercube Sampling approach. This file will output batches of parameter sets to run as .csv files. Each of these .csv files will correspond to a single job distributed to a local High Performance Computing Cluster and the script can be easily modified to change the number of parameters run per job by changing the value of the "per_run" variable.

2. Upload the .csv files for the parameter sets, along with the **"Param_explore.py"** file, to a single directory.

3. Modify the **"Parameter_exploration.slurm"** file so that it sends out the same number of jobs as there are parameter sets, and then run this to distribute the simulations to a local computing cluster.

4. The output "Results" files can be merged together by simply running **"Stitcher.py"** in the same directory as the "Results" files.

5. The "Steady state?" column of the merged results file should be examined to see if any parameter sets fail to reach steady-state.

# Note on modifications

All modifications to how parameters are being calculated should be made in the **"LHS_Sampling_And_Parameter_Export.ipynb"** file. Any and all changes to how fluxes are getting calculated should be made by modifying the "**Param_explore.py**" file. 
