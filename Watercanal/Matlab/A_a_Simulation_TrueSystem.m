close all
clear all
clc

A_b_TrueSystem

Frequency_perturbation = pi/2;  %Frequency of the sinusoidal Perturbation
Amplitude_perturbation = 1;     %Amplitude of the sinusoidal Perturbation
Constant_flow          = 1;     %Magnitud of the constant flow perturbation


dt = 0.01;
Simulation = sim("Sim_Open_Loop");

Animation

% SaveAnimation




