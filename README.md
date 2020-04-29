# Rank epistasis:  A new model for analyzing epistaticinteractions in the absence of quantifiable fitness scores

This repository holds all data from our class project for CSe 845 during MSU's spring 2020 semester. 

Authors:
- Acacia Ackles
- Austin Ferguson
- Connor Grady

This project examines a new, rank-based metric for measuring epistasis in a population. 
In short, the metric exists in four steps:
1. Rank all organisms in population based on performance
2. Mutate all organisms at a particular locus
3. Rank all these mutatnts using the same methods from step 1
4. Compute the edit distance between the two rankings. 

For a thorough dive, we eoncourage you to read the class paper written on the topic here. 

Below is information supplemental to the paper:

## Additional figures
All additional plots are located in the subdirectories of the [./analysis/plots](./analysis/plots) folder. 

## Data analysis scripts
All final analyses are in the root level of the [./analysis/](./analysis) directory. 
Older analyses are still drifiting around, mostly in [./analysis/old](./analysis/old).

## Source code
The only source code that varies from the main branch of MABE is the NKWorld, which is composed of two files located in [./World/NKWorld](./World/NKWorld)

## Replication guide
Thanks to the modularity of [MABE](https://github.com/hintzelab/MABE) and the beaufty of GitHub, replicating the work in this project is simple. 

First we clone the repository:
```
git clone https://github.com/alackles/cse845.git
cd cse845
``` 
Next we need to build MABE:
```
python pythonTools/mbuild.py
```
To run the program, we simply run MABE:
```
./mabe
```
However, we almost always want to run MABE with settings file, so we do that as follows (example settings files included in repo):
```
./mabe -f settings*
```
That's pretty much it!

If you desire to run a large batch of experiments (as we did in the paper), it's easiest to use mq. 
Additional information on mq is given by the help file: 
```python pythonTools/mq.py -h```
mq uses a condition file to specify the runs that will be conducted. Our condition files are included here. 
To run mq with a condition file we run:
```
python pythonTools/mq.py -f {FILENAME}
```
Though to run this you need to either include the flag ```--runLocal``` or ``--runHPCC``. 

That's it! 

If you do run a ton of data, several convenience scripts are included. 
Specifically, ./scrape.R, ./scrape_mutant.R, and ./combine.R which are written to scrape data off MSU's supercomputing cluster with minimal work. 
Then, the data can be passed through the data analyses in [./analysis/](./analysis)
