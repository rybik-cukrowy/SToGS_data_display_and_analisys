This repo is a collection of my nice .root data analisys and display I did during my student internship at IFJ PAN.

For quick overview of contents:
- coincidence.py - plots the coincidence between two tube detectors, where scource energy and their relative orientation is taken from input files. coincidence is where one event appears in two detectors. the code 'splits' the event into two energies and plots a scatter plots of them, where axes are energies deposited in each detector
- timing.py - plots the weighted mean time difference from coincidence events
- rozmycie_czasowe.py - plots scattered broken energy resolution of the hits. Makes simulation data more realistic.
