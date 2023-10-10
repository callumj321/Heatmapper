# Heatmapper
A python module for calculating the contacts to a molecule and then generating an interactive heatmap.

The MDAKits website: https://www.mdanalysis.org/pages/mdakits/

To-do list:

1) Continuous Integration/ Github-hosted runners
2) Test files and folder - short trajectory with 
3) Incorporate RMSD + hydrogen-bonding into heatmapping
4) Add in appendices for the package installation
5) Write to-do list for new features (if necessary)
6) Write jupyter notebook with instructions and notes for new users
7) Change nglview integration to Py3Dmol (https://pypi.org/project/py3Dmol/) for jupyter-notebook integration. iPymol (https://pypi.org/project/ipymol/) might also work.

From meeting 10/10/2023:
1) Eddie and Chris to go through the heatmapper.py code and highlight anything that is not straight forward (CODE REVIEW).
2) Eddie to integrate the RMSD and hydrogen-bonding into heatmapping - need to review how to fit in the code with the new changes to how the colouring/analysis is done, might not work the same with the functions.
