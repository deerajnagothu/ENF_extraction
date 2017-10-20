% STFT anaysis of 50 Hz siganl present in audio signal
%%
clc, close all;
[x, fs] = audioread('Train_Grid_A_A1.wav');  % get the samples of the .wav file

%If the sampling is 44.1kHz or higher then downsampling it by factor of 30
if fs >= 44100
    x = downsample(x,30);
    fs = fs/30;
end

N = length(x);

x = x(:, 1);                                    % get the first channel
x = x(1:300000);         % only first one minute of recording
t = 0:1/fs:(N-1)/fs;
% Plotting the original signal
% figure;
% plot(t,x);   

% plot the fft to check the frequency bins
z = fft(x);
f = (0:N-1)*fs/N;
% figure(2);
% plot(f,abs(z))

figure(10);
z_mags = abs(z);
num_bins = length(z_mags);
% Plotting only the first half of frequencies. That is until nyquist
% frequency
plot((0:1/(num_bins/2 -1):1)*(fs/2), z_mags(1:num_bins/2))


% Spectrogram Parameters 
window=hamming(2048);
noverlap=512;
nfft = 32768;
%nfft=16384;


figure;
spectrogram(x,window,noverlap,nfft,fs,'yaxis')

%% Filtering the signal for the required bandwidth
%Filter Parameters and Plotting Spectrogram
fn = fs/2; %nyquist frequency
Frequency_to_pass = 60;
low_freq = (Frequency_to_pass - 0.5);
high_freq = (Frequency_to_pass + 0.5);

d = designfilt('bandpassiir','FilterOrder',100, ...
   'HalfPowerFrequency1',low_freq,'HalfPowerFrequency2',high_freq, ...
   'SampleRate',fs);


h_designfilt = freqz(d,floor(num_bins/2));
figure(10);
hold on
plot((0:1/(num_bins/2 -1):1), abs(h_designfilt),'g');

y = filter(d,x);

z2 = fft(y);
figure(11);
z2_mags = abs(z2);
num_bins = length(z2_mags);
plot((0:1/(num_bins/2 -1):1), z2_mags(1:num_bins/2))




% figure;
% plot(t,y)
[S,F,T,P]=spectrogram(y,window,noverlap,nfft,fs,'yaxis');
figure;
spectrogram(y,window,noverlap,nfft,fs,'yaxis');

%% Plotting the Time vs Frequency graph. 

[M,N] = size(S);
max_magnitudes = zeros(1,N);
Frequencies = zeros(1,N);
for i = 0:N-1
    max_magnitudes(i+1) = max(abs(S((M*i+1):(M+i*M))));
    location_of_maximum = find(abs(S((M*i+1):(M+i*M))) ==  max_magnitudes(i+1));
    Frequencies(i+1) = F(location_of_maximum,1);
end
% figure
% plot(T,max_magnitudes)

figure;
plot(T,Frequencies)
title('ENF Signal for Linear');

%%  Interpolation of the signal
max_T = max(T);
time = 0:0.5:max_T;
interpolated_frequencies = interp1(T,Frequencies,time,'spline');
figure;
plot(time,interpolated_frequencies)
title('ENF signal Interpolated');

% Examples for filter experimenting
%a = fir1(150,[low_freq high_freq],'bandpass');
%[A,B,C,D] = butter(10,[49 51]/1000,'bandpass');

%sos = ss2sos(A,B,C,D);
% fvt = fvtool(sos,d,'Fs',1000);
% legend(fvt,'butter','designfilt')


% h_butter = freqz(sos,floor(num_bins/2));
% figure(10);
% hold on
% plot((0:1/(num_bins/2 -1):1), abs(h_butter),'r'); 
