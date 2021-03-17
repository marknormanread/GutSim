# GutSim

An agent-based model of the gut microbial community and its repsonse to diet.

## Documentation and help

A comprehensive explanation of how the model works can be [found in this PDF file](doc/gutsim_documentation.pdf).

I launch GutSim from Bash (I use a Mac). 

You can use the `launch.py` script. 
There are several command line options that can be used to supply 
1) the location of the parameters.xml file to use, 
2) the location to which data should be written,
3) the random number seed to use. 

`python launch.py -m <mouse ID here> -p parameters-replicates.xml -d results/example_output_directory -seed 999`

There are more options too. 
You will find a lot of this documented at the bottom of the `main` method at the bottom of `gutsim/experiment.py`.

## Related publications

GutSim was first introduced in:

>> Holmes, A. J., Chew, Y. V., Colakoglu, F., Cliff, J. B., Klaassens, E., Read, M. N., … Simpson, S. J. (2017). Diet-Microbiome Interactions in Health Are Controlled by Intestinal Nitrogen Source Constraints. Cell Metabolism, 25, 140–151. https://doi.org/10.1016/j.cmet.2016.10.021

We are currently working on a follow up paper that focused on the gut microbiome's response to periodic fasting and prebiotics.
