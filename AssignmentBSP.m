%% BSP Assignment
% Last modified: 16-May-19
clear all;
clc;

%% Data Preparation
% The two EEG signals provided (before and during visual stimulation) are
% uploaded to the MATLAB environment
baseline = load('EEG_baseline.mat');
stimulation = load('EEG_stim15Hz.mat');
data_b = baseline.data;
data_s = stimulation.data;

% sampling frequency (equal for baseline and stimulation)
fs = baseline.sampling_rate;
% nyquist frequency
fn = fs/2;
% number of channels (equal for baseline and stimulation)
n_channels = size(baseline.data,1); 

%% Step 1
% 1. Plot in a 5x4 graph the first 5 seconds of the EEG data in the
% baseline condition.

t_max = 5;                  % first five seconds of the recorded data
figure(1)
t = 0:1/fs:(t_max-1/fs);    % horizontal time axis (5s long)

for i = 1:n_channels
    subplot(5,4,i); plot(t,data_b(i,(1:t_max*fs)));
    title(baseline.channels(i).labels);
    xlabel('Time [s]');
    ylabel('Amplitude [mV]');
    xlim([0 5]);
    ylim([-50 50]);
end

%% Step 2
% 2. Plot in a 5x4 graph the Power Spectral Density of the EEG data for the
% baseline condition.

% Anti-aliasing filtering needed? Data already filtered with a bandpass
% 1-45Hz zero phase FIR filter.

% mean removal
for i=1:n_channels
    data_b = data_b - mean(data_b, 2);
end

order = 20; % setting the order of the Yule-Walker estimate
% points of the PSD estimate (zero padding up to the following power of two)
NFFT = 2.^ceil(log2(size(data_b,2)));

% Inizialization of Pxx matrix to save the power spectrum estimation
% for each channel
Pxx_b = zeros(n_channels,NFFT/2+1);

figure(2)
for i = 1:n_channels
    [Pxx_b(i,:), F_b] = pyulear(data_b(i,:),order,NFFT,fs);
    subplot(5,4,i); plot(F_b, 10*log10(Pxx_b(i,:)));
    title(baseline.channels(i).labels);
    ylim([-60 20]);
    xlim([0 fn]);
    xlabel('Frequency [Hz]');
    ylabel('PSD [dB/Hz]');
end

%% Step 3
% 3. Compute the spectral power of each channel in the alpha band which is
% defined in the interval 8-13Hz.

for i=1:n_channels
    % total power of each channel
    p_tot_b(i) = bandpower(Pxx_b(i,:),F_b,'psd');
    % normalized spectral power in the alpha band
    p_alpha_b(i) = bandpower(Pxx_b(i,:),F_b,[8 13],'psd')./p_tot_b(i);
end

figure(3)
topoplot(p_alpha_b',baseline.channels,'maplimits',[min(p_alpha_b) max(p_alpha_b)]);
colorbar();
c = colorbar;
c.Label.String = 'Normalized';
title('Power in the alpha band (baseline condition)');

%% Step 4
% 4. Compute the same as requested in task 3 considering the stimulation
% condition. Compute the change in the alpha power (baseline vs
% stimulation) for each channel.

% To obtain the topoplot of the EEG in the alpha band during stimulation,
% PSD in the stimulation condition must be calculated.

% mean removal
for i=1:n_channels
    data_s = data_s - mean(data_s, 2);
end

% points of the PSD estimate (zero padding up to the following power of two)
NFFT = 2.^ceil(log2(size(data_s,2))); 
% Inizialization of Pxx and F matrix to save the power spectrum estimation
% of each channel
Pxx_s = zeros(n_channels,NFFT/2+1);

% the plot of the PSD is reported in the following task
for i = 1:n_channels
    [Pxx_s(i,:), F_s] = pyulear(data_s(i,:),order,NFFT,fs);
end

% Topoplot under stimulation condition
for i=1:n_channels
    p_tot_s(i) = bandpower(Pxx_s(i,:),F_s,'psd'); % total power of each channel
    p_alpha_s(i) = bandpower(Pxx_s(i,:),F_s,[8 13],'psd')./p_tot_s(i); % normalized spectral power in the alpha band
end

figure(4) % topoplot of the stimulation condition with the same map limits of the baseline one
topoplot(p_alpha_s',stimulation.channels,'maplimits',[min(p_alpha_s) max(p_alpha_s)]);
% the structures stimulation.channels and baseline.channels are equal
colorbar();
c = colorbar;
c.Label.String = 'Normalized';
title('Power in the alpha band (stimulation condition)');

% Compute the change in the alpha power
p_alpha_diff = p_alpha_b - p_alpha_s;

figure(5)
topoplot(p_alpha_diff,stimulation.channels,'maplimits',[min(p_alpha_diff) max(p_alpha_diff)]);
colorbar();
c = colorbar;
c.Label.String = 'Normalized';
title('Power difference in the alpha band');
%% Step 5
% 5. Plot in a 5x4 graph the PSD of the EEG data for the stimulation
% condition.

figure(6)
for i=1:n_channels
    subplot(5,4,i); plot(F_s,10*log10(Pxx_s(i,:)));
    hold on
    xline(8, '-.r');
    hold on
    xline(13, '-.r');
    title(stimulation.channels(i).labels);
    ylim([-60 20]);
    xlim([0 fn]);
    xlabel('Frequency [Hz]');
    ylabel('PSD [dB/Hz]');
    hold off
end

%% Step 6
% 6. Which is the most peculiar change in PSDs when comparing the baselines
% versus the stimulation condition?

% Actually this task does not require any computation.
% However here below PSDs in baseline and stimulation condition are plotted
% overimposed to highlight the differences.

figure(7)
for i=1:n_channels
    subplot(5,4,i);
    plot(F_b,10*log10(Pxx_b(i,:))); 
    hold on
    plot(F_s,10*log10(Pxx_s(i,:)));
    hold on
    xline(8, '-.r');
    hold on
    xline(13, '-.r');
    legend('baseline','stimulation');
    title(stimulation.channels(i).labels);
    ylim([-60 20]);
    xlim([0 fn]);
    xlabel('Frequency [Hz]');
    ylabel('PSD [dB/Hz]');
    hold off
end

% The baseline signals show a peak within the alpha band that does not 
% appear in the stimulation ones; coherent with the decrease of the power 
% in the alpha band already observed.

%%
% 7. After having identified the new frequency originated during
% stimulation, compute the spectral power in the interval [new frequency
% +_1Hz] for each channel. Normalize the derived quantities.

% according to an article this new frequency should be located around 15Hz,
% as the checkerboard switches colors 15 times per second.

% identifying the new frequency
alpha_iso = Pxx_s(:, round(8/fn*(NFFT/2+1)):round(20/fn*(NFFT/2+1)));
[M,I] = max(alpha_iso,[], 2);
I = I + round(8/fn*(NFFT/2+1));
[~, pos] = max(M);
new_freq = F_s(I(pos));

for i=1:n_channels
    % normalized spectral power
    p_new_freq_s(i) = bandpower(Pxx_s(i,:),F_s,[new_freq-1 new_freq+1],'psd')./p_tot_s(i); 
end

figure(8)
topoplot(p_new_freq_s,stimulation.channels,'maplimits',[min(p_new_freq_s) max(p_new_freq_s)]);
colorbar();
c = colorbar;
c.Label.String = 'Normalized';
title('Power in the new stimulation frequency region');

%% Step 8

% Which are the channels in which such new frequency is mostly prominent?
% Can you find in the literature any explanation for this phenomenon? Is
% there any relationship between the power of the new frequency and the
% decrease in alpha power computed in task number 4?

% See report for a detailed answer
