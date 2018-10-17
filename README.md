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



 ### Structure of the Modelscarp_inversion directory
 
```
******** Content of the Modelscarp_inversion folder ******** 


|
└─── Bin							->  folder containing the Modelscarp_inversion executable
|
└─── Data  							->  folder containing the data 

		* data.txt					->  chemical composition of the samples from the fault-planes
		* coll.txt					->  chemical composition of the colluvial wedge
		* sf.txt					->  scaling factors over the time for fast neutrons and muons
											
│   
└─── modelscarp_parameter					
		* param_site.in					->  file containing the parameters of the site and of the inversion
│   
└─── Results						->  folder containing the results files
										
│   
└─── src						->  folder containing the source files
											

```
 

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
    ./configure --prefix=/Users/Jim/Documents/Work/Transdimensional/Modelscarp_RJ-McMC/RJMCMC-4/bin --with-openmpi-include-path /opt/local/include/openmpi-mp --with-openmpi-lib-path=/opt/local/lib/openmpi-mp

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
	Export PKG_CONFIG_PATH = $PKG_CONFIG_PATH:/export/home/RJMCMC/bin/lib/pkgconfig
    ```
   replace "/export/home/RJMCMC/bin/lib/pkgconfig" by your path
    
2. Configure the install file
    ```{r, engine='bash'}
	./configure F77=mpif90 FC=mpif90
	make clean
	make rf_mpi
    ```
    F77= specify the mpi fortran compiler command (here it is "mpif90") 
    FC= specify the mpi #C compiler command (here it is "mpif90") 
    


## How to desing and start an inversion ?

### Searched parameters ?
**Modelscarp Inversion** enables to search the exhumation history that best explain the <sup>36</sup>Cl concentration contained in the bedrock fault-plane. An exhumation history is determined by:

-	The **inheritance history** determined by:
	- the **long-term slip-rate** of the fault that make the samples rise toward the surface before the exhumation of the today observed fault-plane.
	- a potential **quiescence period** of the fault that occurred just prior the exhumation of the today observed fault-plane.
- The **post glacial exhumation history** of the fault-plane that usually includes a part of the fault-plane that has been sampled (the best preserved), and the top-most part of the fault-plane that has not been sampled because the fault-plane is too eroded. This history is parameterized by the **number of exhumation events**, their **ages**, and their **amplitude (slip)**.

<mark>**Modelscarp Inversion** cannot search for the **number of events** since the this parameter has to be fixed for each inversion.</mark> In order to determine the number of events that can be resolved by the <sup>36</sup>Cl approach, **a series of several inversions** has to be achieved with an increasing number of events.

### Inversion project 
The interface enables to design two kinds of inversion project:

- **Single inversion**: a project where a single inversion is operated with a fixed number of event.
- **Multiple inversions**: a project where several inversions are operated with a varying number of event. For instance from 1 to 10 events. 

You can **create a new inversion project** click on `Single inversion`, or  `Multiple inversion`.

<img src="Tutorial/images/Modelscarp_inversion_Single_1.png" width="500">

It also possible to **load an existing project**, in order to modify the parameter of the search by clicking on `Load an inversion project`.

<img src="Tutorial/images/Modelscarp_inversion_existing.png" width="500">


### Create a new project Single inversion
1. Start the assistant by clicking `Single inversion` on the Menu.
2. Fill the parameters value for:
	- the Neighbourhood Algorithm: configuration of the search algorithm.
	- the site parameters: describe the geometry of the site, the density of the bedrock and of the colluvial wedge, and several parameters controlling the <sup>36</sup>Cl production in the samples.
	- the parameters for the scenario of exhumation: the number of events, the parameters you want to search or fix, and the bounds of the search (the pre-exhumation slip-rate, potential quiescence period, age and slip of exhumation events).

	<img src="Tutorial/images/Modelscarp_inversion_Single_2.png" width="500">
3. Import the input files (previously prepared using the excel sheet `data.xls`):
	- chemical data of the samples belonging to the bedrock fault-plane
	- chemical composition of the colluvial wedge
	- neutronic and muonic scaling factors covering the whole duration of the history (including the inheritance, for instance here 300 kyr).
4. Provide a simple name for the inversion project, without symbols!
5. At any time you can save the project to continue later.
6. Generate the inversion files by clicking on `Generate parameter files`.
A new folder for this project is created in the directory `Inversions`. This folder will contains:

```
******** Content of the Inversion folder ******** 

│   
└─── Inversions				->  folder containing all inversion projects 
										
		│   
		└─── TEST					->  Folder of the project 'TEST' 
												
			│   
			└─── bin					->  folder containing the executables 
												"single inversion case"
					- modelscarp_na		->  the inversion program
					- envelop			->  program to determine the
											uncertainties
					- tmp_STOP			-> file to stop an on-going
					 						inversion
					- tmp_misfit		-> file to follow the evolution of
											 the models and misfit during 
											 an inversion
			│   
			└─── data					->  folder containing the data
			│   
			└─── modelscarp_param		->  folder containing the parameter 
											files			
			│   
			└─── results				->  folder containing the results
			
					- na_results.txt	-> result file containing all 
											models tested during the 
											inversion with the associated 
											misfit.
				└─── SITE
					- EQ_sequence.txt	-> final results with the best value 
											found for each parameter and the 
											associated uncertainties.
					- models_in_envelope.txt -> file containing all the 
												models included in the 
												2 sigma uncertainty.	
			│   
			└─── src					->  folder containing the source 
											files
					- makena.macros		-> parameter file containing the 
											compiler commands for the 
											compilation
												
```
### Run the inversion
1. When he inversion files have been generated by clicking on `Generate parameter files`, a new window open. 

	The left part of the window is dedicated to inversion itself: launching the inversion and following the live results (terminal output, and graph of the misfit evolution of tested models over iterations). 

	The right part is dedicated to the plotting of the results after inversion has finished. 
	
	<img src="Tutorial/images/Modelscarp_inversion_Single_3.png" width="600">

2. Start the inversion by clicking on `Start inversion`. This start the program *modelscarp_na*. At anytime, you can stop the inversion by clicking on `Stop inversion`. In that case the inversion will end at the end of the current iteration. This function is particularly interesting if you observe that the misfit has reach a low plateau and is stable, indicating that the best solution has been found. In that case it is not necessary to continue the inversion. 

	<mark>**Important note:** It frequently appears that some tested models have very high misfit value (10E6). Those are models for which the sum of the slip events is larger than the fault-scarp height, thus impossible. For this reason they are directly discarded using a larger misfit value.</mark>
	

3. When inversion has finished, plot the results on the right panel by clicking on `Plot results`. The cumulative slip over the time for each tested model is plotted. The color represent the misfit of the model. The green curve represent the best model. 

	It is also possible to plot the modeled 36Cl profile by clicking on `Plot 36Cl profile`.

	Additional plots can be obtained with `Plot detailed results`.

4. Determine uncertainties associated with the best model by clicking `Determine uncertainties`. This starts the program *envelop* that can take a certain time. The evolution can be followed in the terminal output. When uncertainties are determined, it is possible to plot the cumulative slip over the time of the models included within the analytical uncertainty. 

5. Results from an inversion can also be plotted latter by clicking `Plot results` of the main menu.

### Create a new project of Multiple inversion
In that case, the project will be composed of several inversion with a number of events varying between a minimum and maximum value that you determine.

1. Start the assistant by clicking `Multiple inversion` on the Menu.

2. Similarly to the single inversion, fill the different parameters. 
	
	<img src="Tutorial/images/Modelscarp_inversion_mult_1.png" width="600"> 
	
	For the exhumation scenario, provide the minimum and the maximum number of events that you want to test. For instance if you choose a minimum of 1 and a maximum of 3 events, the interface will create 3 different inversions, one with a single exhumation event (called *1ev*), a second with 2 events (called *2ev*), and a third with 3 events...

	In that case, 3 folders are created in the inversion directory : *1ev*, *2ev*, and *3ev*, each containing the executables for the inversion. See the example below:
	
	
	```
******** Content of the Inversion folder ********    
└─── Inversions		->  folder containing all inversion projects 
										
		│   
		└─── TEST		->  Folder of the multi-inversion project 'TEST' 
			│   
			└─── 1ev				-> folder containing the executables    
				└─── bin				for the inversion using 1 event
				└─── data	   
				└─── modelscarp_param   
				└─── results	   
				└─── src					
			│   
			└─── 2ev				-> folder containing the executables    
				└─── bin				for the inversion using 2 events
				└─── data	   
				└─── modelscarp_param   
				└─── results	   
				└─── src	
			│   
			└─── 3ev				-> folder containing the executables    
				└─── bin				for the inversion using 3 events
				└─── data	   
				└─── modelscarp_param   
				└─── results	   
				└─── src	
				
	```
3. `Generate parameter files for the inversions` create those inversion folders and start the assistant to run the inversion and plot the live results. 

4. Because the inversion can take a long time on a single computer, it is recommended to run the multiple inversions out of Matlab. The Matlab graphical interface will enable you to follow the results. 

	To run an inversion:
	
   1. open a terminal (out of Matlab)

   2. Go to the inversion directory. For instance for the 1 event inversion:

		```	
		cd Modelscarp_V2/Inversions/TEST/1ev/bin
		```
		
   3. Start the inversion with the following command:
		
		```	
 		./modelscarp
 		```	
 		
5. Because the inversion can take a long time on a single computer, it is recommended to run the multiple inversions out of Matlab. The Matlab graphical interface will enable you to follow the results by clicking on `Live results`, or on `Plot results (on-going or finished inversion)` on the initial menu.

	<img src="Tutorial/images/Modelscarp_inversion_mult_2.png" width="300"><img src="Tutorial/images/Modelscarp_inversion_mult_3.png" width="200">

6. This interface enables you to follow the progress of each inversion that you have previously started, by selecting which inversion is followed on the left panel.

	<img src="Tutorial/images/Modelscarp_inversion_mult_4.png" width="500">
	
	-> When an inversion has finished, the results can be plotted on the right panel (cumulative slip over the time of each tested model, the color of the curve varies as function of the misfit of the model)

    <img src="Tutorial/images/Modelscarp_inversion_mult_5.png" width="500">

	-> Detailed results can be plotted showing for each tested models, the age and the slip of the events and the associated misfit (RMSw).

	<img src="Tutorial/images/Modelscarp_inversion_mult_6.png" width="500">
	
7. Once an inversion has finished, you need to run the program ```envelop```  in the terminal (not in Matlab) to find all the models included within the analytical uncertainty and determine the uncertainties on the exhumation scenario. 

	-> The cumulative slip of all those models are figured by the light green curves, and the best model in dark green.

	<img src="Tutorial/images/Modelscarp_inversion_mult_7.png" width="500">
		
	-> It is possible to plot the <sup>36</sup>Cl profile of the best model (red dots), compared with the data (black dots), as shown by the following figure:	

	<img src="Tutorial/images/Modelscarp_inversion_mult_8.png" width="300">
	
8. When all inversions have finished, it is possible to plot the RMSw of the best models has function of the number of events (as shown in the following figure):

	<img src="Tutorial/images/Modelscarp_inversion_mult_9.png" width="300">

9. It is also possible to plot the cumulative slip of the final scenario inferred from each inversion with a varying number of events:	

	<img src="Tutorial/images/Modelscarp_inversion_mult_10.png" width="400">
	
## How to run modelscarp on a cluster ?

It is better to run inversion on a cluster because usually a large number of iterations are required to converge and well explore the parameter space. 

### Design the inversion project
1. Similarly to the previous example, you need to create a `multiple inversion` project. In that case, **you had to specify that the inversion will be operated on a cluster** by clicking on `cluster`. Then generate the parameter files for the inversions.

	<img src="Tutorial/images/Modelscarp_inversion_cluster.png" width="400">

2. **Copy the folder of the inversions on the cluster.**
3. **The executables of each inversion have to be compiled on the cluster** by running those commands in the terminal of your cluster, located in the directory of the inversion project.
	
	For instance:

	1. Connection from my computer to the cluster with the terminal:

		```	
		ssh yourcluster.com -l login
		```
	
	2. Go to the directory of the inversion project called "TEST":

		```	
		cd your_path_to_modelscarp/Inversions/TEST
		```
	
	3. Give the right to execute the script:

		```	
		chmod u+rwx compile_all.sh
		```
	
	4. Execute the script of compilation:

		```	
		./compile_all.sh
		```
		wait for one or two minutes...
		
4. **Start each inversion using the appropriated command of your cluster.** For instance for the 1 event inversion, the script may contain those commands:

	```	
	% go to the binaries folder of the 1 event inversion
	cd your_path_to_modelscarp/Inversions/TEST/1ev/bin
	
	% run the inversion program
	./modelscarp_na	
	
	% run the program to determine the uncertainties
	./envelop
	```
4. **At any time you can follow the evolution of an inversion** by downloading the inversion folder on your computer, and using the Matlab interface `Plot results (on-going or finished inversion)` (in  the first menu). 
	In that case, the file `tmp_misfit`in the `bin` directory will be used to follow the evolution of the misfit over the iteration. 
	
	<img src="Tutorial/images/Modelscarp_inversion_mult_5.png" width="400">
	
	___
