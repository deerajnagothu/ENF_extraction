% STFT anaysis of 50 Hz siganl present in audio signal
%%
clc, close all;
[x, fs] = audioread('Train_Grid_A_P1.wav');  % get the samples of the .wav file

prime_frequency = 60;
enf = 1; % 0 for power signals 1 for audio signals
%If the sampling is 44.1kHz or higher then downsampling it by factor of 30
if fs >= 44100
    x = downsample(x,30);
    fs = fs/30;
end



x = x(:, 1);                                    % get the first channel
x = x(1:60000);         % only first one minute of recording
N = length(x);
t = 0:1/fs:(N-1)/fs;
% Plotting the original signal
% figure;
% plot(t,x);   

% plot the fft to check the frequency bins
z = fft(x);
f = (0:N-1)*fs/N;
% figure(2);
% plot(f,abs(z))


%figure(10);
z_mags = abs(z);
num_bins = length(z_mags);
% Plotting only the first half of frequencies. That is until nyquist
% frequency

%plot((0:1/(num_bins/2 -1):1), z_mags(1:num_bins/2))


% Spectrogram Parameters 
window=hamming(2048);
noverlap=1600;
nfft = 32768;
%nfft=16384;


figure;
spectrogram(x,window,noverlap,nfft,fs,'yaxis')
% Collecting each frequency separately
sig_prime1 = zeros(1,N);
sig_prime2 = zeros(1,N);
sig_prime3 = zeros(1,N);
sig_prime4 = zeros(1,N);
sig_prime5 = zeros(1,N);
sig_prime6 = zeros(1,N);

upper_limit_frequency = 6* prime_frequency;
%% Filtering the signal for the required bandwidth
%Filter Parameters and Plotting Spectrogram
for fq = prime_frequency:prime_frequency:upper_limit_frequency
    fn = fs/2; %nyquist frequency
    Frequency_to_pass = fq;
    low_freq = (Frequency_to_pass - 0.5);
    high_freq = (Frequency_to_pass + 0.5);

    d = designfilt('bandpassiir','FilterOrder',100, ...
       'HalfPowerFrequency1',low_freq,'HalfPowerFrequency2',high_freq, ...
       'SampleRate',fs);


    h_designfilt = freqz(d,floor(num_bins/2));
    %figure(10);
    %hold on
    %plot((0:1/(num_bins/2 -1):1), abs(h_designfilt),'g');

    y = filter(d,x);

    z2 = fft(y);
    %figure(11);
    z2_mags = abs(z2);
    num_bins = length(z2_mags);
    %plot((0:1/(num_bins/2 -1):1), z2_mags(1:num_bins/2))




    % figure;
    % plot(t,y)
    [S,F,T,P]=spectrogram(y,window,noverlap,nfft,fs,'yaxis');
    %figure;
    %spectrogram(y,window,noverlap,nfft,fs,'yaxis');

    % Plotting the Time vs Frequency graph. 

    [M,N] = size(S);
    max_magnitudes = zeros(1,N);
    Frequencies = zeros(1,N);
    
    for i = 0:N-1
        max_magnitudes(i+1) = max(abs(S((M*i+1):(M+i*M)))); 
        location_of_maximum = find(abs(S((M*i+1):(M+i*M))) ==  max_magnitudes(i+1));
        
        Frequencies(i+1) = F(location_of_maximum,1);
        if fq == prime_frequency
            sig_prime1(i+1) = Frequencies(i+1);
        end
        if fq == 2*prime_frequency
            Frequencies(i+1) = Frequencies(i+1) - prime_frequency;
            sig_prime2(i+1) = Frequencies(i+1);
        end
        if fq == 3*prime_frequency
            Frequencies(i+1) = Frequencies(i+1) - 2*prime_frequency;
            sig_prime3(i+1) = Frequencies(i+1);
        end
        if fq == 4*prime_frequency
            Frequencies(i+1) = Frequencies(i+1) - 3*prime_frequency;
            sig_prime4(i+1) = Frequencies(i+1);
        end
        if fq == 5*prime_frequency
            Frequencies(i+1) = Frequencies(i+1) - 4*prime_frequency;
            sig_prime5(i+1) = Frequencies(i+1);
        end
        if fq == 6*prime_frequency
            Frequencies(i+1) = Frequencies(i+1) - 5*prime_frequency;
            sig_prime6(i+1) = Frequencies(i+1);
        end
    end
    % figure
    % plot(T,max_magnitudes)
    color = ['k','r','g','b','y','m'];
    j = fq/prime_frequency;
    col = color(j);
    figure(5);
    hold on
    plot(T,Frequencies,col)
    title('ENF Signal for Linear');

    %  Interpolation of the signal
    max_T = max(T);
    time = 0:0.5:max_T;
    interpolated_frequencies = interp1(T,Frequencies,time,'spline');
    figure(6);
    hold on
    plot(time,interpolated_frequencies,col)
    title('ENF signal Interpolated');
end

figure(5)
hold off
legend('60Hz','120Hz','180Hz','240Hz','300Hz','360Hz')


figure(6)
hold off
legend('60Hz','120Hz','180Hz','240Hz','300Hz','360Hz')

% Removing zero's from the end of the signal

sig_prime1 = sig_prime1(1:N);
sig_prime2 = sig_prime2(1:N);
sig_prime3 = sig_prime3(1:N);
sig_prime4 = sig_prime4(1:N);
sig_prime5 = sig_prime5(1:N);
sig_prime6 = sig_prime6(1:N);

figure(11);
plot(T,sig_prime1);
hold on
plot(T,sig_prime3);
hold on
plot(T,sig_prime5);
ylabel('60Hz,180Hz & 300 Hz Signal');
xlabel('Time');
hold off
legend('60Hz','180Hz','300Hz');

figure(12);
plot(T,sig_prime2);
hold on
plot(T,sig_prime4);
hold on
plot(T,sig_prime6);
ylabel('120Hz,240Hz & 360 Hz Signal');
xlabel('Time');
hold off
legend('120Hz','240Hz','360Hz');

%% Cross-Correlation between the signals
if enf == 0
    s1 = sig_prime1;
    s2 = sig_prime3;
    s3 = sig_prime5;
else
    s1 = sig_prime2;
    s2 = sig_prime4;
    s3 = sig_prime6;
end

%Testing by adding a little delay
%s2 = [zeros(1,50) s2(1:end-50)];


%computing the cross-correlation between the three pairs of signals and
%normalize them such that the maximum value is one.

[C21,lag21] = xcorr(s2,s1);
C21 = C21/max(C21);

[C31,lag31] = xcorr(s3,s1);
C31 = C31/max(C31);

[C32,lag32] = xcorr(s3,s2);
C32 = C32/max(C32);

% location of max values indicate time leads or lags
 
[M21,I21] = max(C21);
t21 = lag21(I21);

[M31,I31] = max(C31);
t31 = lag31(I31);

[M32,I32] = max(C32);
t32 = lag31(I32);

% Plotting the correlation with the lags
figure(7);
subplot(3,1,1)
plot(lag21,C21,[t21 t21],[-0.5 1],'r:')
text(t21+100,0.5,['Lag: ' int2str(t21)])
ylabel('C_{21}')
axis tight
title('Cross-Correlations')

subplot(3,1,2)
plot(lag31,C31,[t31 t31],[-0.5 1],'r:')
text(t31+100,0.5,['Lag: ' int2str(t31)])
ylabel('C_{31}')
axis tight

subplot(3,1,3)
plot(lag32,C32,[t32 t32],[-0.5 1],'r:')
text(t32+100,0.5,['Lag: ' int2str(t32)])
ylabel('C_{32}')
axis tight
xlabel('Samples')


