# nmrR1Z_sim
--------------------------------------------------------------------------------------
Package for calculating carbon-hydrogen (CH) bond R1Z relaxation rates from all-atom MD simulations
 
Protocol based on Doktorova, Khelashvili, Ashkar and Brown, 2022, Biophysical Journal, 122:984-1002  
https://doi.org/10.1016/j.bpj.2022.12.007

If you use the results from the code in a publication, please cite the article above.

--------------------------------------------------------------------------------------
To run the calculation, need the following in a single directory:

- _traj_center.dcd_: trajectory centered so that mean z position of all terminal methyl groups in the bilayer is at z=0. Example script for centering the trajectory is included in the additional_scripts/ directory. This is not critical but helps with keeping everything consistent.

- _struct.psf_: psf file for the trajectory

- _box2.txt_: file with simulation box dimensions. Example script for calculating the box dimensions is included in the additional_scripts/ directory.

- nmrR1 directory included in this package: a separate nmrR1 directory will be needed for each lipid for which the CH bond dynamics need to be analyzed. For example, in mixtures the directory can be renamed to _nmrR1_resname_ where 'resname' is the name of the lipid


--------------------------------------------------------------------------------------
In the nmrR1 directory, update the _ch_vectors.txt_ file by editing:

- the name of the lipid to be analyzed (resname in VMD)
- the names of the lipid's carbon and hydrogen atoms for which CH bond dynamics will be calculated; each row should have 3 atom names, 1 carbon followed by 2 hydrogen atoms. The names of the atoms will be used to make atom selections. For double bonds where a carbon has only one hydrogen atom, write the name of the hydrogen atom twice. For example, for the oleoyl chain of POPC, using CHARMM36 notation, the file will include:

C28 H8R H8S  
C29 H9R H9R  
C210 H10R H10R  
C211 H11R H11S  

--------------------------------------------------------------------------------------
To calculate the relaxation rates:

1. Get the direction of all CH vectors in every frame of the trajectory with the following:

```
./run_analysis > log.out &
```

This may take a while depending on the size of the trajectory; start process in the background 

2. Get the CH bond autocorrelation functions and order parameters from the output files from Step 1. To do this, open MATLAB inside the directory with all output files from Step 1 and run the _calculate_ACF_Scd.m_ script. MATLAB can be started from the terminal without a GUI by typing:

```
[path_to_matlab]/matlab -nodesktop
```

The script can then be executed by simply typing its name:

```
calculate_ACF_Scd
```

3. Calculate the relaxation rates from the ACFs with the _calculate_R1.m_ script from within MATLAB (open the file first and modify the user-specified parameters if necessary):

```
calculate_R1
```

-------------------------------------------
To visualize and further analyze the data:

- Plot R1 vs order parameter squared for different subsets of carbon atoms. Example script is provided in nmrR1/results_dtt/_plot_squarelaw.m_

- Fit the data to a line and use its slope to calculate the bending rigidity of the bilayer. Example script is provided in nmrR1/results_dtt/_get_kappa_from_squarelaw.m_

- To estimate the full bilayer thickness from the number density profile of the bilayer, calculate the number density of the whole bilayer (example in additional_scripts/_calculate_ndp.tcl_) and find the distance between the slabs where the NDP falls below 5% of the max value. For example, in MATLAB that can be done with the following commands:

```
ndp = load('ndp.dat');  
maxleft = max(ndp(1:floor(length(ndp(:,1))/2),2));  
maxright = max(ndp(floor(length(ndp(:,1))/2):end,2));  
maxpeak = mean([maxleft maxright]);  
ind = find(ndp(:,2)<0.05*maxpeak);  
indleft = ind(ind<length(ndp(:,1))/2);  
indright = ind(ind>length(ndp(:,1))/2);  
full_bilayer_thickness = ndp(indright(1),1)-ndp(indleft(end),1)
```
