# LaggedRecovery
Code for VAST model and additional analyses presented in Robertson et al. (In Review.)


The analyses for this manuscript occurred in four parts which are represented by separate scripts.

Part 1 involved spatial bottom water temperature kriging, the code for this is located in full_temp_krig_process.R

Part 2 involved running a VAST model for american plaice and yellowtail flounder, this script is located in vast_model_run.R

Part 3 involved running a non-linear mixed effects model for the residuals from the yellowtail flounder VAST, this script is located in exponential_fxn.R and exponential_fxn_RE.cpp

Finally, part 4 can some aspects of part 3 can be found in the markdown file, post_vast_analyses.md

This github also contains the data files that should be necessary to run every aspect of the analyses except the initial VAST due to data availability restrictions imposed by Fisheries and Oceans Canada. These data can be requested by following the instructions in the Data Availability section of the manuscript. Upon request I (mdrobertson24@gmail.com) can provide the code used to prepare that data for the VAST model.
