
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all 
close all

% Set Working Directory

cd 'C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\Rust Simulation'

% clear the raw data

clear_data

% estimate the tables

make_tables

% generate the montecarlo simulation

montecarlo
