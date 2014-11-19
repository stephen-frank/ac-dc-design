The `Batch` directory contains templates for batch execution of the monolith or
NGBD algorithm for the mixed AC-DC electricity distribution system design model.
The scripts are intended to be modified by the user as needed for particular
runs of the model.

* `batch-run-monolith.gms` executes the monolith with various global solvers.
* `batch-run-algorithm.gms` executes the NGBD algorithm.
* `solve-acdc-model.gms` and `check-acdc-model.gms` are helper scripts used by
  the batch scripts. These must be in the same directory as the batch scripts.
  
To perform a batch solve for the model:

1. If desired, copy the batch scripts to a new directory for the run.
2. Copy the appropriate data file(s) from the `Data` directory to the local
   directory.
3. Modify `batch-run-algorithm.gms` and/or `batch-run-monolith.gms` as needed.
   a. Specify the GAMS command
   b. Specify the path to your GAMS license file
   c. Perform any other required configuration changes.  
      (See descriptive comments in the _Setup_ section of the batch script.)
4. Execute `batch-run-algorithm.gms` and/or `batch-run-monolith.gms`.

The batch scripts output the solved model to `[Filename]_solved.gdx`, in which
`[Filename]` is as specified in the script.

See also the README and the edited batch scripts in the `Case_Studies` directory
for examples.
