File descriptions:  

LHS_Sampling_And_Parameter_Export.ipynb: This is a Jupyter notebook which conducts parameter sampling and generates .csv files containing, in each row of each file, a unique set of parameters to run in a simulation.  

Param_explore.py: This is a Python script that iterates through all of the parameter sets contained in a single .csv file, runs the corresponding simulations, and outputs the results (along with a flag of whether or not the simulation reached steady state).  

Parameter_exploration.slurm: This UNIX file distributes the individual param_explore.py jobs to the computing cluster. 

Stitcher.py: This Python file collects all of the individual results files together and turns them into a single .csv.  

Using these files 

First, by modifying the loop depth and individual parameter value lists in the Jupyter notebook, you can make your parameter list .csv files. Note that the per_run variable determines how many parameter sets (i.e. how many rows) to have in each file. Based on runs I did on my local system, I think 400 is a reasonable number. If the total number of parameter sets is not evenly divisible by 400, the last .csv file with just have the remainder.  
Modify parameter_exploration.slurm file (ideally, make a copy and rename) so that it sends out the appropriate number of jobs. Note that the .csv files exported by the parameter sheet exporter are 1-indexed. You will modify the “#SBATCH --array=1-100” line so that the last number equals the number of .csv files you are going to iterate through.  
Modify stitcher.py; just change the second number in the line for i in range(2,5): to be whatever the final results file number will be.  
Once you have generated all of your individual .csv parameter lists, you need to upload these files, the param_explore.py file, and the parameter_exploration.slurm file to the cluster. I recommend putting all of these files in a dedicated folder in your scratch drive so as to keep everything contained and so you aren’t using up your main directory’s size limit.  
When everything is in place, login to the computing cluster and navigate to the folder where you have uploaded your files.  
Run:  
sbatch parameter_exploration.slurm   
Once results are generated, you can navigate back to the scratch folder. Now we are going to run stitcher.py to collect everything into one big file. Run the following two commands: 
module load Conda/3 
python stitcher.py 
With the big dataset in hand, you can now move that final file onto your local drive for downstream analysis.  
The first thing I would check after loading into your analysis platform of choice would be the “steady state?” column – if any “Fail” values are present, this means that one or more simulations failed to reach steady-state. 
