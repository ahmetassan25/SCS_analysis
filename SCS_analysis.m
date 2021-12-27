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

%% Filter data
fl = 100;
fu = 500;
fs = data.streams.EMGw.fs;

[b,a] = butter(2, [(2*fl)/fs (2*fu)/fs],'bandpass');
final_biceps = filtfilt(b,a,biceps);
final_ecr = filtfilt(b,a,ecr);
final_fcr = filtfilt(b,a,fcr);

%% extract the stimulation/response window
% cd 'G:\.shortcut-targets-by-id\1DLmzFkqnKdtnmtAAuam9z99HmOFaf-OA\scs_mapping\scgs002\scgs002_2021-11-30_scs-implant\ephys\scgs002_2021-11-30_01_scsmap'
stim_time = data.scalars.eS1p.ts;
direct = dir('*.mat');
stim_table = load (direct.name);
stim_num = size(stim_table.T);
num_of_stim = stim_num(1); % Opens table to check the total number of stim

stim_intensities = unique(stim_table.T.pulse_amplitude);
count_int = histc(stim_table.T.pulse_amplitude,stim_intensities);

snap_biceps = {};
snap_ecr = {};
snap_fcr = {};
window_len = 0.05;
for i=1:num_of_stim; 
    
    snap_biceps = [snap_biceps [final_biceps((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];
    snap_ecr = [snap_ecr [final_ecr((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];
    snap_fcr = [snap_fcr [final_fcr((stim_time(i)-window_len)*fs : (stim_time(i)+window_len)*fs)]];

end

%% plot
tt = numel(snap_biceps{1})/fs;
time = -tt/2:1/fs:(tt/2-(1/fs));

%% plot all lines

yy = zeros(numel(snap_biceps{1}),1);

figure
for i=1:numel(snap_biceps);
    plot(time, snap_biceps{i}+yy')
    hold on
    yy = yy+.01;
end

%% plot only the same amplitudes
select_amp = 200;
test_amp = snap_fcr(stim_table.T.pulse_amplitude==select_amp);
amp_ord = 0;
for i=1:numel(test_amp);
    plot(time,test_amp{i}+amp_ord);
    hold on
    amp_ord = amp_ord+0.01
end

%% plot same stim configurations

%organize/combine the channels
ch1 = stim_table.T.channel1;
ch2 = stim_table.T.channel2;

all_ch1 = {};
all_ch2 = {};
ch_combined_all = {};
for i=1:length(ch1)
    ch1_new = '';
    ch2_new = '';
    cm = ',';
    for ii=1:4;
        if ii == 4;
            cm = '';        
        end
        ch1_new = strcat(ch1_new, num2str(ch1(i,ii)),cm);  
        ch2_new = strcat(ch2_new, num2str(ch2(i,ii)),cm);  
    end
    ch_combined = strcat (ch1_new,',',ch2_new);
    all_ch1{end+1} = ch1_new;
    all_ch2{end+1} = ch2_new;    
    ch_combined_all{end+1} = ch_combined;     
end

all_ch1 = all_ch1';
all_ch2 = all_ch2';
ch_combined_all = ch_combined_all';
ch_unique = unique(ch_combined_all);

table_test = stim_table.T;
table_test.ch_all = ch_combined_all;

% select stimulation pairs and plot vs intensity
 %%
ch_uniques = unique(table_test.ch_all);
fprintf('there are %d different stimulation pair',numel(ch_uniques))
unique_check = [];
unique_stim_protocol = 12; % select an unique stim protocol.
for i=1:length(ch_combined_all)
    unique_check = [unique_check strcmp(table_test.ch_all{i},ch_uniques{unique_stim_protocol})];
end

indx_stim = find(unique_check==1); %find the location of selected stim pair
if strcmp(table_test.ch_all{1}, ch_uniques{unique_stim_protocol})==0 % check if the baselin is included
    indx_stim = [1 indx_stim]; %include the first stim with 0 mA as a control
end

amp_ord = 0;
figure;
subplot(1,3,1)
gap_plots = 0.0001;
for i=indx_stim;
    subplot(1,3,1)
    plot(time,snap_biceps{i}+amp_ord);hold on;
    subplot(1,3,2)
    plot(time,snap_ecr{i}+amp_ord);hold on;
    subplot(1,3,3)
    plot(time,snap_fcr{i}*.1+amp_ord);
    amp_ord = amp_ord+gap_plots;
    hold on;
end


subplot(1,3,1);box off;title('Biceps');ylabel('Amplitude (uA)');ylim([-2e-4 18e-4 ])
yticklabels(num2cell(unique(table_test.pulse_amplitude)))
yticks([0:gap_plots:(gap_plots*(numel(indx_stim)-1))])
subplot(1,3,2);box off;title('ECR');yticklabels(num2cell([]));xlabel('Time (ms)');ylim([-2e-4 25e-4 ])
subplot(1,3,3);box off;title('FCR');yticklabels(num2cell([]));ylim([-2e-4 25e-4 ])

%% AUC 
auc_biceps = [];
auc_ecr = [];
auc_fcr = [];
% for i=1:numel(snap_biceps);
%     auc_biceps = [auc_biceps sum(abs(snap_biceps{i}(1245:1428)))];% 15*fs/1000 (15ms time window)
%     auc_ecr = [auc_ecr sum(abs(snap_ecr{i}(1245:1428)))];
%     auc_fcr = [auc_fcr sum(abs(snap_fcr{i}(1245:1428)))];
% end

for i=1:numel(snap_biceps);
    auc_biceps = [auc_biceps sum(abs(snap_biceps{i}(1245/2:1428/2)))];% 15*fs/1000 (15ms time window)
    auc_ecr = [auc_ecr sum(abs(snap_ecr{i}(1245/2:1428/2)))];
    auc_fcr = [auc_fcr sum(abs(snap_fcr{i}(1245/2:1428/2)))];
end

amps = unique(table_test.pulse_amplitude);
% figure
% plot(amps , (auc_biceps(indx_stim)/auc_biceps(1))*100-100,'.-');hold on;
% plot(amps,(auc_ecr(indx_stim)/auc_ecr(1))*100-100,'.-');hold on;
% plot(amps,(auc_fcr(indx_stim)/auc_fcr(1))*100-100,'.-');
% % norm_auc_fcr = auc_fcr(indx_stim)/max(auc_fcr(indx_stim))

figure
plot(amps , (auc_biceps(indx_stim)/auc_biceps(1))*100,'.-');hold on;
plot(amps,(auc_ecr(indx_stim)/auc_ecr(1))*100,'.-');hold on;
plot(amps,(auc_fcr(indx_stim)/auc_fcr(1))*100,'.-'); box off
% % norm_auc_fcr = auc_fcr(indx_stim)/max(auc_fcr(indx_stim))

hold on;
x_labels = linspace(amps(1),amps(end),numel(amps));
xticklabels(num2cell(amps))
xticks(amps)
%% Threshold
% find(norm_auc(1)*3<norm_auc)
norm_fcr = (auc_fcr(indx_stim)/auc_fcr(1))*100-100;
threshold_fcr = find(norm_fcr>norm_fcr(1)*2)
norm_ecr = (auc_ecr(indx_stim)/auc_ecr(1))*100-100;
threshold_ecr = find(norm_ecr>norm_ecr(1)*2)
norm_biceps = (auc_biceps(indx_stim)/auc_biceps(1))*100-100;
threshold_biceps = find(norm_biceps>norm_biceps(1)*2)
all_norms = [norm_fcr; norm_ecr]
max_auc = max(cat(2,norm_fcr, norm_ecr,norm_biceps))
figure;plot(ones([0:max_auc],1)*amps(4),[0:max_auc])

% hold on; plot(ones(numel([0:max_auc]))*amps(threshold_ecr(1)),[0:max_auc])
% plot(ones(numel([0:max_auc]))*amps(threshold_fcr(1)),[0:max_auc])
% plot(ones(numel([0:max_auc]))*amps(threshold_biceps(1)),[0:max_auc])




