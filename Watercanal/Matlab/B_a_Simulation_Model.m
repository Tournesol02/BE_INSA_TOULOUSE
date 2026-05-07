close all
clear all
clc

B_b_Model
Frequency_perturbation = pi/2;
Amplitude_perturbation = 1;
Constant_flow          = 1;

dt = 0.01;
Simulation = sim("Sim_Open_Loop");

Animation
