# Welcome to SISYPHOS
<font color='red'>**THIS MODULE IS IN DEVELOPMENT, ATM YOU CANNOT STOP A RUNNING PROCEDURE OTHER THAN KILLING OLEX2**</font>
<br>
<br>


This tool is meant to provide means of serially refine multiple models against multiple datasets of the same structure.
In other words, it solves Sisyphos dilemma for you, repeating the same task for the same structure over and over again. 
It can be used for finding out which crystal was the best out of many (looking at synchrotron users!), or for benchmarking theoretical calculations. 

<br>

NoSpherA2 (* Kleemiss et al., Chem. Sci., 2021, 12, 1675
&nbsp; URL[https://pubs.rsc.org/en/content/articlehtml/2021/sc/d0sc05526c,PAPER] ) is implemented and can be used in this procedure.

<br>
<br>

Anomalous dispersion refinement (* Meurer et al., IUCrJ, 2022, 9 
&nbsp; URL[https://journals.iucr.org/m/issues/2022/05/00/lt5050/,PAPER]) can be also used, especially when in the vicinity of absorption edges (see Anomalous Dispersion in Info tab).

# Inputs
## Set working dir

Here you can show the program the location of your measured data. SISYPHOS will look for any .hkl that is in any subfolder here, except for output folders of SISYPHOS.


## Set solution ins

Here you can give a starting model in form of an .ins file. Further filetypes such as .cif and .res will be implemented in the future. The model WILL be refined for every dataset, so this is really just the starting point.

# Adjustment
## Adjustment / eV

This setting is only for synchrotron users who know an offset between the enery / wavelength in their headers and the exact energy / wavelength of their experiment. Leave this empty if you measured on an in-house difffractometer or you don't care about minute differences in the wavelength.

# Resolution
## Resolution / Angstroms

Here you can apply a general cutoff in d-spacing for all of your datasets.


# Elements
Decide which Elements the Anomalous Dispersion values should be refined for.


# disp
If you have more than one position for one of the Elements given, you can choose if they should have the same DISP value (constrained) or an individual one for each site.

