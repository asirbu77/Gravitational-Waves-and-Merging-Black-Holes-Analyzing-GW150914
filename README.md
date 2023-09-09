# Gravitational-Waves-and-Merging-Black-Holes-Analyzing-GW150914

## Project description
This project was conducted as a course requirement for Astronomy 4602 at Western University.

## How to use the code
1 - RUN SignalProcessing.py, it will output two csv files into your working directory, ligoDat.csv and
    NRDat.csv (These csv files are already available in the repository, but they are also updated by the code file)

2 - RUN Frequency_v_Mass_Sim.py to obtain frequency-mass relation of black hole mergers.

3 - RUN MassFinder_4points.py and MassFinder_12points.py to determine the masses of the black holes and see plots.

4 - RUN TheoreticalMerger.py to model the inspiral of the black holes based on Newtonian mechanics

5 - RUN EquilibriumFunction.py to investigate validity of compactness ratio of 1.7

## Scientific report
The report can be found in the repository as Gravitational-Waves-and-Merging-Black-Holes-Analyzing-GW150914.pdf

## Notes:

Frequency_v_Mass_Sim.py, MassFinder_4points.py, MassFinder_12points.py, RUN TheoreticalMerger.py, and EquilibriumFunction.py
were all coded by Aidan Sirbu and Katie Brown.

SignalProcessing.py was repurposed from https://www.gw-openscience.org/GW150914data/GW150914_tutorial.html with slight
modifications made to it. readligo.py is a supplementary file used by SignalProcessing.py. 
