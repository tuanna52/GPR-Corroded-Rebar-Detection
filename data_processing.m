%This file pre-processes the data imported from “Steel_In_Sandstone_Raw_Data_Only.mat” file.
%Then applies inverse model to derive physical parameters of each layer inside a structure

load Steel_In_Sandstone_Raw_Data_Only.mat;

E0 = 1; %Amplitude of transmitted signal
fs = 1700; %Sampling frequency of the system (Hz);
F0 = 24e9; %Start frequency of the chirp (Hz)
Bc = 1.5e9; %Bandwidth of the chirp (Hz)
Tc = 300e-3; %Sweep time of the chirp (s)

%Extract and process signal from the raw data
x = TestLiveX_See_ppt(2:end, 1);
x = x(:)';
x = normalize(x,'center','mean');
x = x/10000;

t = (0:length(x)-1)/fs;
frequency = fs*(0:(length(t)/2-1))/length(t);

%Low-pass filter to remove frequency under 20 Hz
E_if = lowpass(x, 20, fs);

%Plot data after pre-processing
plot(t, E_if);
xlabel("time (s)");
ylabel("Amplitude (V/m)");
title("Data after pre-processing")

%Fourier transform
ps2 = abs(fft(E_if))/length(t);

figure;
plot(frequency(1:100),ps2(1:100));
xlabel("Frequency (Hz)");
ylabel("P1(f)");
title("Fourier transform")

%Applying inverse model to derive physical parameters of each layer
figure;
para_estimate = inverse_model(E_if, E0, fs, F0, Bc, Tc);