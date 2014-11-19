The `Test` directory contains a set of test scripts for the mixed AC-DC
electricity distribution system design formulation and solution algorithm. The
tests build sequentially. To run the full suite of tests:

1. Copy the appropriate data file from the `Data` directory to the `Test`
   directory. (The default configuration requires the file `threebus.gdx`.)
2. Ensure the relative path to the AC-DC model and algorithm files (in the
   `Formulation` folder) is set correctly for your local configuration.
3. Execute the seven test scripts in the following order:
   a. `test-sets.gms`
   b. `test-data.gms`
   c. `test-big-M.gms`
   d. `test-bound-tightening.gms`
   e. `test-model.gms`
   f. `test-solve.gms`
   g. `test-algorithm.gms`
   
The tests allow a number of configuration options; these are documented at the
top of each test script. In addition, some test scripts have special
instructions (also documented in the script). For example, a complete test of
the sets requires running `test-sets.gms` twice with slightly different options.
