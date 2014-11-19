The `Case_Studies` directory contains preconfigured case studies which match
the test cases presented in the article:

> S. M. Frank and S. Rebennack, "Optimal Design of Mixed AC-DC Distribution
> Systems for Commercial Buildings: A Nonconvex Generalized Benders
> Decomposition Approach," *European Journal of Operational Research*, 2014,
> accepted for publication.

To prepare for the case studies, unzip the entire electronic appendix to your
local file system, preserving relative paths.

To execute a particular case study:

1. Copy the appropriate data file from the `Data` directory to the case study
   directory. (To save space, these files are pre-copied.)
2. Edit the batch run script to match your local GAMS configuration:
   a. Specify the GAMS command
   b. Specify the path to your GAMS license file
   c. Perform any other required configuration changes.  
      (See descriptive comments in the _Setup_ section of the batch script.)
3. Execute the batch script using GAMS:
   a. `batch-run-monolith.gms` for the monolith
   b. `batch-run-algorithm.gms` for the NGBD algorithm
   
The batch scripts output the solved model to `[Filename]_solved.gdx`, in which
`[Filename]` is as specified in the script.
