# Modelscarp inversion
   Modelscarp Inversion  Copyright (C) 2017  TESSON J. and BENEDETTI L. 2017

[![DOI](https://zenodo.org/badge/99244458.svg)](https://zenodo.org/badge/latestdoi/99244458)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed for research purposes in the hope that
   it will be useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
   the GNU General Public License for more details.

   Use it on your own risk.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For any question, report bugs, propose improvements..., please contact:
        - Tesson J. at jim.tesson@gmail.com, or
        - Benedetti L. at benedetti@cerege.fr

Constraining the past seismic activity and the slip-rates of faults over several millennials is crucial for seismic hazard assessment. Chlorine 36 (36Cl) in situ produced cosmogenic nuclide is increasingly used to retrieve past earthquakes histories on seismically exhumed limestone normal fault-scarps. Following Schlagenhauf et al., [2010] modelling approach, we present a new methodology to retrieve the exhumation history based on a Bayesian transdimensional inversion of the 36Cl data. This procedure uses the reversible jump Markov chains Monte-Carlo algorithm (RJ-MCMC, Green [1995]) which enables 1-exploring the parameter space (number of events, age and slip of the events), 2-finding the more probable scenarios, and 3-precisely quantifying the associated uncertainties. Through a series of synthetic tests, the algorithm revealed a great capacity to constrain event slips and ages in a short computational time (several days) with a precision that can reach 0.1 ky and 0.5 m for the age and slip of exhumation event, respectively. In addition, our study show that the amount of 36Cl accumulated when the sampled fault-plane was still buried under the colluvial wedge, prior  its exhumation, might represents up to 35 % of the total 36Cl. This contribution can be accurately determined with a depth profile, reducing uncertainty on the exhumation scenario.


Please cite the use of this code using the following DOI:

[![DOI](https://zenodo.org/badge/99244458.svg)](https://zenodo.org/badge/latestdoi/99244458)


## Getting Started

The inversion of the 36Cl data is based on a modified version of the RJMcMC algorithm provided by Gallagher et al. (2011) (http://www.iearth.org.au/codes/rj-MCMC/). It is thus required to first install the RJMCMC library, and then install the inversion routine.

### Prerequisites

To use Modelscarp Inversion, you will need:

- an Unix system (OSX, Linux...)
- an MPI library (e.g. openMPI)
- a C compiler (e.g. gcc, https://gcc.gnu.org)
- a Fortran compiler (e.g. gfortran that is included in gcc)
- Makefile (https://www.gnu.org/software/make/)

 

### How to install the library RJ-McMC

The RJ-McMC library is provided in the folder RJMCMC. The installation is operated using your terminal.

1.  Specify the path for the mpi library in the terminal
	for instance: 
    ```{r, engine='bash'}
    cd RJMCMC
    PATH=$PATH:/export/apps/mpich2/bin
    ```
    replace "/export/apps/mpich2/bin" by the absolute path of your mpi library.

2. Configure the install file
	for instance: 
    ```{r, engine='bash'}
    ./configure --prefix=/Users/Jim/RJMCMC/bin --with-openmpi-include-path /opt/local/include/openmpi-mp --with-openmpi-lib-path=/opt/local/lib/openmpi-mp

    ```
	Options:
	
--prefix : specify the path of the bin folder where the RJ-McMC library will be installed.

--with-openmpi-lib-path= specify the path of the “lib” folder of your mpi library.

--with-openmpi-include-path= specify the path of “include” folder of your the mpi library.


3. Install the RJ-McMC library : 
    ```{r, engine='bash'}
	make clean
	make
	make install
    ```

###  Installation of the program Modelscarp Inversion

1.Specify the "pkgconfig" path localized in the directory of RJMCMC (bin/lib/pkconfig) 
 
    ```{r, engine='bash'}
	cd Modelscarp_Inversion
	Export PKG_CONFIG_PATH = $PKG_CONFIG_PATH:/Users/Jim/RJMCMC/bin/lib/pkgconfig
    ```
   replace "/Users/Jim/RJMCMC/bin/lib/pkgconfig" by your path
    
2. Configure the install file
    ```{r, engine='bash'}
	./configure F77=mpif90 FC=mpif90
	make clean
	make rf_mpi
    ```
    F77= specify the mpi fortran compiler command (here it is "mpif90") 
    
    FC= specify the mpi #C compiler command (here it is "mpif90") 
    
 ### Structure of the Modelscarp_inversion directory
 
```
******** Content of the Modelscarp_inversion folder ******** 


|
└─── Bin							->  folder containing the Modelscarp_inversion executable
|
└─── Data  							->  folder containing the data 

		* data.txt					->  chemical composition of the samples from the fault-planes
		* coll.txt					->  chemical composition of the colluvial wedge
		* EL.txt					->  geomagnetic scaling factors over the time for fast neutrons and muons
											
│   
└─── modelscarp_parameter					
		* param_site.in					->  file containing the parameters of the site and of the inversion
│   
└─── Results						->  folder containing the results files
										
│   
└─── src						->  folder containing the source files
											

```

## How to start an inversion ?

### Searched parameters ?
**Modelscarp Inversion** enables to search the exhumation scenarios that best explain the <sup>36</sup>Cl concentration contained in the bedrock fault-plane. An exhumation history is determined by:

-	The **inheritance history** determined by:
	- the **long-term slip-rate** of the fault that make the samples rise toward the surface before the exhumation of the today observed fault-plane.
	- a potential **quiescence period** of the fault that occurred just prior the exhumation of the today observed fault-plane.
- The **post glacial exhumation history** of the fault-plane that usually includes a part of the fault-plane that has been sampled (the best preserved), and the top-most part of the fault-plane that has not been sampled because the fault-plane is too eroded. This history is parameterized by the **number of exhumation events**, their **ages**, and their **amplitude (slip)**.


### Inverse your dataset
1. Prepare the data files in the “data” folder, (following the excel sheet `Format_your_data.xls`):
	- data.txt : chemical data of the samples belonging to the bedrock fault-plane
	- coll.txt : chemical composition of the colluvial wedge
	- EL.txt : neutronic and muonic scaling factors covering the whole duration of the history (including the inheritance, e.g. 300 000 yr + 20 000 yr = 320 000 yr. The EL file must cover at least 320 000 yr).

2. Edit the parameter file “modelscarp_param/param_site.in:
    ```{r, engine='bash'}
	####################################################
	#  		PARAMETERS OF SITE
	#			
	####################################################
	####################################################
	#       Input data files parameters
	####################################################
	#### colluvium data file: coll.txt
	1 : number of line
	62 : number of column  ( ** do not change !)
	#### fault-scarp data file: data.txt
	101 : number of line
	66 : number of column ( ** do not change !)
	#### scaling factors file: sf.txt
	5001 : number of line
	4 : number of column  ( ** do not change !)
	####################################################
	#       Site Parameters
	####################################################
	800 : Total height of the post-glacial scarp (cm)
	200 : Maximum sample depth (cm below the colluvial wedge surface)
	0.0 : Erosion rate (mm/yr)
	#### Fault-scarp geometry
	20 : alpha (colluvial wedge surface angle)
	50 : beta (fault scarp surface angle)
	30 : gamma (upper peri-glacial surface angle)
	#### Density
	2.00 : colluvium
	2.70 : rock
	####################################################
	#       Elementary production rates
	####################################################
	42.2 : spallation rate in Ca
	2.303e-06 : 	lambda_36 Radioactive decay constant for 36Cl (a-1)
	208 : 	True attenuation length for fast neutron (g.cm-2)
	####################################################
	#       Inversion Parameters
	####################################################
	#### Slip-rate prior exhumation
	y : search ? (y/n) if no, it is fixed with the minimum value
	0.0 : min slip-rate (mm/yr)
	5.0 : max slip-rate (mm/yr)
	0.5 : Std dev. of slip-rate value change
	#### Length of the long-term history prior the post-glacial exhumation
	300000 : duration (yr)
	#### Quiescence period
	y : Include a quiescence period prior the exhumation of the scarp ? (y/n)
	50000 : max quiescence period length if it is searched (yr)
	2000.0 : Std dev. of scarp top age (yr)
	#### Number of events
	3 : Minimum number of events
	20 : Maximum number of events
	#### Event ages
	0.0 : Min age
	20000.0 : Max age
	#### Algorithm search parameters
	20.0 : Std dev. of move change (pd) (cm)
	300.0 : Std dev. of age value changes (yr)
	300.0 : Standart dev. of birth/death events (yr)
	1200000 : number of iteration total
	101 :	seed
	983 : seed mult

    ```
    

3. Prior the execution of the Modelscarp_Inversion program, specify in the terminal the path for the library of the RJ-McMC in the bin folder (bin/lib):
    ```{r, engine='bash'}
	cd Modelscarp_inversion
	export LD_LIBRARY_PATH = $LD_LIBRARY_PATH:/Users/Jim/RJMCMC/bin/lib
    ```
4. Start the Modelscarp_Inversion program:
    ```{r, engine='bash'}
	./Modelscarp_inversion
    ```
5. On-going results for each chain of the inversion are placed in the folder 'results'.
    ```{r, engine='bash'}
	ChainNumber_results.txt
    ```
	e.g: for the chain 2, the result file is named: “2_results.txt”
	
	Each line of the file is a step of the chain (= an accepted model). On a line you will find:

    ```{r, engine='bash'}
	Number of events, Slip-rate (mm/yr) Quiescence period length (kyr, if searched), event ages, event slips, rmsw.
    ```
    
	e.g: the following model is composed of 2 events at 8291 and 9047 yr with a slip of 997 and 602 cm, a peri-glacial slip-rate of 3.52 mm/yr and no quiescence period. The RMSw of the model is 333.78.
	
    ```{r, engine='bash'}
	2, 3.52, 8291 9047, 997.253505 602.746495, 333.78
    ```

6. At the end of the inversion, the results of each chain are avalaible in the folder 'results'. 

