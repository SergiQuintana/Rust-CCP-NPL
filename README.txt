README
Author: Sergi Quintana Garcia
-------------------------------

This package generates the results of Rust (1987) “Optimal Replacement of GMC Bus Engines” using full solution methods, CCPs, and Nested Pseudo-Likelihood from Aguirregabiria and Mira (2002). It also performs a Monte Carlo exercise to compare the performance of the different methods with varying sample sizes.

------------------------
Instructions for Users:
------------------------

To produce all outputs, run main.m. All outputs will be stored in the Output folder.

Codes:

--Clear_data.m: Clears the Rust (1987) data.
--Make_tables.m: Estimates the different methods (full solution, CCPs, Aguirregabiria and Mira) using the Rust (1987) data.
--Montecarlo.m: Performs the Monte Carlo exercise.

Outputs:

--Transitions.tex: Estimates of the transition probabilities.
--Thetas_full.tex: Estimates of the thetas using full solution methods.
--Allmethods.tex: Estimates of the thetas using all three methods.
--Montecarlo_table.tex: Estimates of the thetas from the Monte Carlo simulation exercise.

Other Codes in the Folder:

--Estimation_rust.m: Estimates the Rust model with Rust data. Run after clear_data.m.
--Estimation_rustccp.m: Estimates the Rust model using CCPs with Rust data. Run after clear_data.m.
--Estimation_aguirre.m: Estimates the Rust model using the Aguirregabiria and Mira approach with Rust data. Run after clear_data.m.
--Simulation_data.m: Simulates and stores data from the Rust model. Run before Estimation_full.
--Estimation_full.m: Estimates the Rust model using the simulated data. Run after Simulation_data.m.
--Estimation_ccps.m: Estimates the Rust model using CCPs with simulated data. Run after Simulation_data.m.
--Estimation_aguirregabira.m: Estimates the Rust model using the Aguirregabiria and Mira method with simulated data. Run after Simulation_data.m.
