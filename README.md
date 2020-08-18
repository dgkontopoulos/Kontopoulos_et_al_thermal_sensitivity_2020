
This repository contains code that reproduces the analyses of the following study:

>Dimitrios - Georgios Kontopoulos, Thomas P. Smith, Timothy G. Barraclough, and Samraat Pawar: **[Adaptive evolution shapes the present-day distribution of the thermal sensitivity of population growth rate](https://doi.org/10.1371/journal.pbio.3000894)**. PLOS Biology. 2020.

<br>

---

#### Dependencies

The scripts are written in **shell**, **Python 2.7.6**, **Perl 5.22.1**, and **R 3.4.4**. Besides the interpreters, the following modules or packages are needed:

* **For Python 2.7.6:**
    * bigfloat, version 0.2.1
	* joblib, version 0.11
	* lmfit, version 0.8.1
	* numpy, version 1.14.1
	* rpy2, version 2.8.2
	* scipy, version 0.14.0<br><br>
* **For Perl 5.22.1:**
 	* Statistics::R, version 0.34<br><br>
* **For R 3.4.4:**
	* ade4, version 1.7-13
	* ape, version 5.2
	* BAMMtools, version 2.1.6
	* boot, version 1.3-20
	* coda, version 0.19-2
	* cowplot, version 0.9.3
	* geiger, version 2.0.6.4
	* MCMCglmm, version 2.26
	* motmot.2.0, version 1.1.2
 	* msm, version 1.6.6
 	* phytools, version 0.6-60
 	* R.devices, version 2.16.1
 	* Rphylopars, version 0.2.12
 	* spptest, version 0.4

Other required software:

* BayesTraitsV3, version 3.0.2
* levolution 
* stabletraits and stabletraitssum

---
 
#### Execution

Scripts have to be executed from within `Code`. For R scripts, in particular, some of them produce R objects that need to be manually examined. Thus, R scripts should generally be run from within R, using the `source` command. Data files needed for execution (see https://doi.org/10.6084/m9.figshare.12816140.v1) have to be placed in a `Data` directory, outside `Code`. Outputs will be written in a `Results` directory, outside `Code`.
