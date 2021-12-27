test
clear all
close all
clc

%% Load Data
cd 'G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs002\scgs002_2021-12-2_Day(2)_scs_stim(0-200)\ephys'

addpath(genpath(' C:\Users\asa2248\Desktop\mapping\TDTMatlabSDK (1)\TDTSDK\TDTbin2mat')); %Add with subfolders
addpath(genpath('G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs002\scgs002_2021-12-2_Day(2)_scs_stim(0-200)'))
data = TDTbin2mat('G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs002\scgs002_2021-12-2_Day(2)_scs_stim(0-200)\ephys')

%data = TDTbin2mat('C:\Users\asa2248\Desktop\mapping\scgs001_2021-11-09-01-cxs_mu-20211109T193500Z-001\scgs001_2021-11-09-01-cxs_mu')
emg_data = data.streams.EMGw.data; 
biceps = double(emg_data(1,:));
ecr = double(emg_data(2,:));
fcr = double(emg_data(3,:));
