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
8) Get the code to work for residues, atoms, indices.

From meeting 10/10/2023:
1) Eddie and Chris to go through the heatmapper.py code and highlight anything that is not straight forward (CODE REVIEW).
2) Eddie to integrate the RMSD and hydrogen-bonding into heatmapping - need to review how to fit in the code with the new changes to how the colouring/analysis is done, might not work the same with the functions.
3) Callum to review/update doc strings and write high-level/broad text (few sentences) in the README to discuss what the code does.
4) Chris to look into the test files etc. and review if they need restructuring.

From meeting 01/11/2023:
1) Eddie and Callum to decide on what systems we should be producing to use as test cases. Melittin for dimerising, lipid contacts; ubiquitin for protein-protein(?) interactions [lots of binding domains] 
2) Eddie to chase Micaela for a second meeting
3) Chris to write up what he needs for testing
4) Survey users for requirements spec - document names: requirement specification, Moscow analysis (must have, could have, will not have)
