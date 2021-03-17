# GutSim
An agent-based model of the gut microbial community and it's repsonse to diet

## Documentation and help

A comprehensive explanation of how the model works can be (found in this PDF file)[doc/gutsim_documentation.pdf].

I launch GutSim from Bash (I use a Mac). 

You can use the `launch.py` script. 
There are several command line options that can be used to supply 
1) the location of the parameters.xml file to use, 
2) the location to which data should be written,
3) the random number seed to use. 

`python launch.py -m <mouse ID here> -p parameters-replicates.xml -d results/example_output_directory -seed 999`

There are more options too. 
You will find a lot of this documented at the bottom of the `main` method at the bottom of `gutsim/experiment.py`.

