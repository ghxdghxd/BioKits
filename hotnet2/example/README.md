# HotNet2 Examples

We provide example config files for creating influence matrices, running `simpleRun.py`, as well as
generating heat scores, delta selection, and running Hotnet2 with/without significance testing.
Below we demonstrate how to run HotNet2 using these examples. Please consult the config files in
`example/configs` to see the parameters used for each example.

## File organization
All examples are located in the `example/` directory. The example directory by default includes
test data for running HotNet2 using mutation data or a heat file, as well as a test network.

## Create influence matrices
First, you must create the influence matrix and directory of permtued influence matrices used by
HotNet2 in these examples. To do so, run:

    python makeRequiredPPRFiles.py @example/configs/influence_matrix.config
    
## A Simple Run
To run the simpleRun example, execute:

    python runHotNet2.py @example/configs/simple.config
    
This will run delta selection, HotNet2, and statistical significance testing using the heat file
`example/example.heat` and store the output in `example/output/simple`. Because delta is selected
using a permutation test, the output may vary slightly from run to run. You can see visualizations
of the subnetworks in `example/output/simple/viz`, e.g.

![](http://f.cl.ly/items/1V0i2S2U3G0m003l0H0N/Screen%20Shot%202014-01-10%20at%204.47.14%20PM.png).

## Advanced usage examples
To run an example using mutation data, first create the directory `example/output/mutation`. Then:

1. Generate the heat file from the mutation data.
       
        python generateHeat.py @example/configs/heat.config

2. Perform delta selection.

        python bin/findThreshold.py @example/configs/delta.config
        
   (Note that `findThreshold.py` will output a file that contains a distribution of deltas. In
   practice, you would use this distribution to choose the delta for running HotNet2 (e.g. by
   taking the median), and then update the `delta` parameter of
   `example/configs/significance.config` and `example/configs/run.config`)

3. Run HotNet2 (without significance testing).

        python bin/findComponents.py @example/configs/run.config
        
4. Run HotNet2 (with significance testing).

        python bin/findComponents.py @example/configs/significance.config

5. Make results website.

		python bin/makeResultsWebsite.py @example/configs/website.config
