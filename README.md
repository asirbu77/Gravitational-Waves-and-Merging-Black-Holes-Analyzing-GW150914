# Gravitational-Waves-and-Merging-Black-Holes-Analyzing-GW150914

## Abstract
The purpose of this investigation is to explore the properties of gravitational waves and the evolution of black hole mergers. This is performed in the context of GW150914 using a classical model with minimal general relativistic corrections. The gravitational wave event was detected by the LIGO on September 14, 2015 at 09:50:45 UTC. The source of GW150914 has been established in previous work to be a binary black hole merger. The determination of the individual masses is based on Newtonian and Keplerian dynamics. This investigation explores the effectiveness of classical physics in predicting the behaviour and characteristic masses of the black holes. Firstly, the relationship between gravitational wave frequency and black holes mass is determined. The gravitational wave frequency at merger is used to determine a total mass of $68.36M_\odot$. A linear regression of $f_{GW}^{-8/3}$ vs time is used to find a chirp mass of $27.81M_\odot$. Using varying ranges of frequency sample points, the masses of the two black holes were determined to lie in the ranges of $42.4M_\odot \text{ to } 45.35M_\odot$ and $23.02M_\odot \text{ to } 25.96M_\odot$. These classically determined masses are in moderate agreement with cited values deduced from general relativistic equations. Furthermore, Keplerian orbital equations as well as Newtonian mechanics are used to model the time-evolution of the system's orbital frequency, velocity and separation.

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
