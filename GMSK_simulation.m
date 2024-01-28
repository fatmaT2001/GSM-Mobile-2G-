
clear; clc; close all;

%% processing the image and converting it into a binary stream

% Read the image file
originalImage = imread('photo.jpg');  % Replace with your image path
% Convert the image to grayscale (if it's not already)
grayscaleImage = rgb2gray(originalImage);

% Display the original grayscale image
figure(1);
imshow(grayscaleImage);
title('Original Image');

% Convert the grayscale image to a binary stream
% First, convert the image to double for processing
doubleImage = double(grayscaleImage);
% Convert each pixel value to binary
binaryImage = de2bi(doubleImage);
% Reshape the binary matrix to a binary stream (vector)
binarySignal = reshape(binaryImage', [], 1);

%% creating the object required for GMSK modulation 
% Initialize GMSK modulator and demodulator
gmskModulator = comm.GMSKModulator('BitInput', true, 'PulseLength', 1, 'SamplesPerSymbol', 1);
gmskDemodulator = comm.GMSKDemodulator('BitOutput', true, 'PulseLength', 1, 'SamplesPerSymbol', 1);

%% modulating the signal and initialize needed parameters 
% Modulate the binary signal using GMSK
modulatedSignal = gmskModulator(binarySignal);

% Define simulation parameters
velocity = 100 * 1e3 / 3600; % Travelling speed in m/s
carrierFrequency = 920e6; %Carrier frequency in Hz
maxDopplerShift = velocity * carrierFrequency / physconst('lightspeed'); % Maximum Doppler shift
sampleRate = 270.833e3; % Symbol rate in symbols per second
SNR_dB = 1:1:40; % SNR range in dB

%% create the object required for error calculation
% Initialize error rate calculator
errorRateCalculator = comm.ErrorRate('ReceiveDelay', gmskDemodulator.TracebackDepth);

%% AWG channel
% Initialize an array to store BER for each SNR value
BER_AWGN = [];

% Preparing to display received images at different SNR levels on the same figure
numImagesToDisplay = sum(mod(1:length(SNR_dB), 5) == 0); % Counting how many images we will display
figReceivedImages = figure; % Create a figure for all received images

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Process the modulated signal through the AWGN channel
    noisySignal = awgnChannel(modulatedSignal);

    % Demodulate the noisy signal
    demodulatedSignal = gmskDemodulator(noisySignal);

    % Calculate the error statistics
    errorStats = errorRateCalculator(binarySignal, demodulatedSignal);
    BER_AWGN = [BER_AWGN errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));

        % Determine the subplot position for the current image
        subplotIndex = 1 + (i - 5) / 5;
        subplot(2, ceil(numImagesToDisplay / 2), subplotIndex);
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% Adjusting the layout and title of the figure with received images
sgtitle(figReceivedImages, 'Received Images at Different SNR Levels (AWG)');

% Plot the BER vs SNR for the AWGN channel
figure(2);
semilogy(SNR_dB, BER_AWGN);
title('BER vs SNR plot (AWGN Channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Rayleigh channel flat fast fading
% Initialize an array to store BER for the Rayleigh channel
BER_ray = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rayleigh channel object for flat fast fading
    rayleighChannel = comm.RayleighChannel(...
        'SampleRate', sampleRate, 'PathDelays', [0 1e-20], 'AveragePathGains', [2 3],...
        'NormalizePathGains', true, 'MaximumDopplerShift', maxDopplerShift, 'RandomStream', 'mt19937ar with seed',...
        'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rayleigh channel
    rayleighSignal = rayleighChannel(awgnChannel(modulatedSignal));

    % Demodulate the Rayleigh channel processed signal
    demodulatedSignal = gmskDemodulator(rayleighSignal);

    % Calculate the error statistics
    errorStats =errorRateCalculator(binarySignal, demodulatedSignal);
    BER_ray = [BER_ray errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(3)
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Flat Fast Fading Rayleigh channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the Rayleigh flat fast fading channel
figure(4);
semilogy(SNR_dB, BER_ray);

title('BER vs SNR plot (Flat Fast Fading Rayleigh channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();
%%  Rayleigh channel frequency-selective fast fading
% Initialize an array to store BER for the frequency-selective fast fading Rayleigh channel
BER_ray2 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rayleigh channel object for frequency-selective fast fading
    rayleighChannel = comm.RayleighChannel(...
        'SampleRate', sampleRate, 'PathDelays', [10.0 8.0]*1e-6, 'AveragePathGains', [2 3],...
        'NormalizePathGains', true, 'MaximumDopplerShift', maxDopplerShift, 'RandomStream', 'mt19937ar with seed', ...
        'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rayleigh channel
    rayleighSignal = rayleighChannel(awgnChannel(modulatedSignal));

    % Demodulate the Rayleigh channel processed signal
    demodulatedSignal = gmskDemodulator(rayleighSignal);

    % Calculate the error statistics
    errorStats = errorRateCalculator(binarySignal, demodulatedSignal);
    BER_ray2 = [BER_ray2 errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(5)
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Frequency Selective Fast Fading Rayleigh channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the frequency-selective fast fading Rayleigh channel
figure(6);
semilogy(SNR_dB, BER_ray2);
title('BER vs SNR plot (Frequency Selective Fast Fading Rayleigh channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% %% Rayleigh channel flat slow fading
% Initialize an array to store BER for the Rayleigh flat slow fading channel
BER_ray3 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rayleigh channel object for flat slow fading
    rayleighChannel = comm.RayleighChannel(...
        'SampleRate', sampleRate, 'PathDelays', [0 1e-20], 'AveragePathGains', [2 3],...
        'NormalizePathGains', true, 'MaximumDopplerShift', 0, 'RandomStream', 'mt19937ar with seed', ...
        'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rayleigh channel
    rayleighSignal = rayleighChannel(awgnChannel(modulatedSignal));

    % Demodulate the Rayleigh channel processed signal
    demodulatedSignal = gmskDemodulator(rayleighSignal);

    % Calculate the error statistics
    errorStats = errorRateCalculator(binarySignal, demodulatedSignal);
    BER_ray3 = [BER_ray3 errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(7)
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Flat Slow Fading Rayleigh Channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the Rayleigh flat slow fading channel
figure(8);
semilogy(SNR_dB, BER_ray3);
title('BER vs SNR plot (Flat Slow Fading Rayleigh channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Rayleigh channel frequency-selective slow fading
BER_ray4 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rayleigh channel object for frequency-selective slow fading
    rayleighChannel = comm.RayleighChannel(...
        'SampleRate', sampleRate, 'PathDelays', [10.0 8.0]*1e-6, 'AveragePathGains', [2 3],...
        'NormalizePathGains', true, 'MaximumDopplerShift', 0, 'RandomStream', 'mt19937ar with seed',...
        'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rayleigh channel
    rayleighSignal = rayleighChannel(awgnChannel(modulatedSignal));

    % Demodulate the Rayleigh channel processed signal
    demodulatedSignal = gmskDemodulator(rayleighSignal);

    % Calculate the error statistics
    errorStats = errorRateCalculator(binarySignal, demodulatedSignal);
    BER_ray4 = [BER_ray4 errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(9)
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Frequency Selective Slow Fading Rayleigh channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the frequency-selective slow fading Rayleigh channel
figure(10);
semilogy(SNR_dB, BER_ray4);
title('BER vs SNR plot (Frequency Selective Slow Fading Rayleigh channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Rician Channel flat fast fading
BER_ric = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannel = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rician channel object for flat fast fading
    ricianChannel = comm.RicianChannel(...
        'SampleRate', sampleRate, 'PathDelays', [0.0 1.0]*1e-20, 'AveragePathGains', [2 3],...
        'KFactor', 2.8, 'NormalizePathGains', true, 'DirectPathDopplerShift', 5.0, ...
        'DirectPathInitialPhase', 0.5, 'MaximumDopplerShift', maxDopplerShift, 'DopplerSpectrum', doppler('Bell', 8),...
        'RandomStream', 'mt19937ar with seed', 'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rician channel
    ricianSignal = ricianChannel(awgnChannel(modulatedSignal));

    % Demodulate the Rician channel processed signal
    demodulatedSignal = gmskDemodulator(ricianSignal);

    % Calculate the error statistics
    errorStats = errorRateCalculator(binarySignal, demodulatedSignal);
    BER_ric = [BER_ric errorStats(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(11)
        binaryMatrix = reshape(demodulatedSignal, 8, [])'; 
        decimalMatrix = bi2de(binaryMatrix);
        reconstructedImage = reshape(decimalMatrix, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Flat Fast Fading Rician channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImage));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the flat fast fading Rician channel
figure(12);
semilogy(SNR_dB, BER_ric);
title('BER vs SNR plot (Flat Fast Fading Rician channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Rician Channel flat slow fading
BER_ric2 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannelFlatSlow = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rician channel object for flat slow fading
    ricianChannelFlatSlow = comm.RicianChannel(...
        'SampleRate', sampleRate, 'PathDelays', [0.0 1.0]*1e-20, 'AveragePathGains', [2 3],...
        'KFactor', 2.8, 'NormalizePathGains', true, 'MaximumDopplerShift', 0, 'DopplerSpectrum', doppler('Bell', 8),...
        'RandomStream', 'mt19937ar with seed', 'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rician channel
    ricianSignalFlatSlow = ricianChannelFlatSlow(awgnChannelFlatSlow(modulatedSignal));

    % Demodulate the Rician channel processed signal
    demodulatedSignalFlatSlow = gmskDemodulator(ricianSignalFlatSlow);

    % Calculate the error statistics
    errorStatsFlatSlow = errorRateCalculator(binarySignal, demodulatedSignalFlatSlow);
    BER_ric2 = [BER_ric2 errorStatsFlatSlow(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(13)
        binaryMatrixFlatSlow = reshape(demodulatedSignalFlatSlow, 8, [])'; 
        decimalMatrixFlatSlow = bi2de(binaryMatrixFlatSlow);
        reconstructedImageFlatSlow = reshape(decimalMatrixFlatSlow, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Flat Slow Fading Rician Channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImageFlatSlow));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the flat slow fading Rician channel
figure(14);
semilogy(SNR_dB, BER_ric2);
title('BER vs SNR plot (Flat Slow Fading Rician channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();
%% Rician channel (Frequency-selective fast fading)
BER_ric3 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannelFreqSelectiveFast = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rician channel object for frequency-selective fast fading
    ricianChannelFreqSelectiveFast = comm.RicianChannel(...
        'SampleRate', sampleRate, 'PathDelays', [10.0 8.0]*1e-6, 'AveragePathGains', [2 3],...
        'KFactor', 2.8, 'NormalizePathGains', true, 'DirectPathDopplerShift', 5.0, 'DirectPathInitialPhase', 0.5, ...
        'MaximumDopplerShift', maxDopplerShift, 'DopplerSpectrum', doppler('Bell', 8), 'RandomStream', 'mt19937ar with seed', ...
        'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rician channel
    ricianSignalFreqSelectiveFast = ricianChannelFreqSelectiveFast(awgnChannelFreqSelectiveFast(modulatedSignal));

    % Demodulate the Rician channel processed signal
    demodulatedSignalFreqSelectiveFast = gmskDemodulator(ricianSignalFreqSelectiveFast);

    % Calculate the error statistics
    errorStatsFreqSelectiveFast = errorRateCalculator(binarySignal, demodulatedSignalFreqSelectiveFast);
    BER_ric3 = [BER_ric3 errorStatsFreqSelectiveFast(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(15)
        binaryMatrixFreqSelectiveFast = reshape(demodulatedSignalFreqSelectiveFast, 8, [])'; 
        decimalMatrixFreqSelectiveFast = bi2de(binaryMatrixFreqSelectiveFast);
        reconstructedImageFreqSelectiveFast = reshape(decimalMatrixFreqSelectiveFast, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Frequency Selective Fast Fading Rician Channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImageFreqSelectiveFast));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the frequency-selective fast fading Rician channel
figure(16);
semilogy(SNR_dB, BER_ric3);
title('BER vs SNR plot (Frequency Selective Fast Fading Rician channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Rician Channel frequency-selective slow fading
BER_ric4 = [];

% Loop through each SNR value
for i = 1:length(SNR_dB)
    % Create an AWGN channel object with the specified SNR
    awgnChannelFreqSelectiveSlow = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)', 'SNR', SNR_dB(i));

    % Create a Rician channel object for frequency-selective slow fading
    ricianChannelFreqSelectiveSlow = comm.RicianChannel(...
        'SampleRate', sampleRate, 'PathDelays', [10.0 8.0]*1e-6, 'AveragePathGains', [2 3],...
        'KFactor', 2.8, 'NormalizePathGains', true, 'MaximumDopplerShift', 0, 'DopplerSpectrum', doppler('Bell', 8),...
        'RandomStream', 'mt19937ar with seed', 'Seed', 22, 'PathGainsOutputPort', true);

    % Process the modulated signal through the Rician channel
    ricianSignalFreqSelectiveSlow = ricianChannelFreqSelectiveSlow(awgnChannelFreqSelectiveSlow(modulatedSignal));

    % Demodulate the Rician channel processed signal
    demodulatedSignalFreqSelectiveSlow = gmskDemodulator(ricianSignalFreqSelectiveSlow);

    % Calculate the error statistics
    errorStatsFreqSelectiveSlow = errorRateCalculator(binarySignal, demodulatedSignalFreqSelectiveSlow);
    BER_ric4 = [BER_ric4 errorStatsFreqSelectiveSlow(1)]; % Append the BER for this SNR

    % Reconstruct and display the received image at certain SNR levels
    if mod(i, 5) == 0
        figure(17)
        binaryMatrixFreqSelectiveSlow = reshape(demodulatedSignalFreqSelectiveSlow, 8, [])'; 
        decimalMatrixFreqSelectiveSlow = bi2de(binaryMatrixFreqSelectiveSlow);
        reconstructedImageFreqSelectiveSlow = reshape(decimalMatrixFreqSelectiveSlow, size(grayscaleImage));
        subplot(2, ceil(length(SNR_dB) / 10), i / 5); % Plotting in the same figure
        sgtitle('Frequency Selective Slow Fading Rician Channel Received Images')
        str = sprintf('SNR =%d dB', i); 
        imshow(uint8(reconstructedImageFreqSelectiveSlow));
        title(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    end
end

% BER plot for the frequency-selective slow fading Rician channel
figure(18);
semilogy(SNR_dB, BER_ric4);
title('BER vs SNR plot (Frequency Selective Slow Fading Rician channel)');
xlabel('SNR in dB');
ylabel('BER');
ylim([0 0.5]);
grid on;

% Release system objects to free resources
errorRateCalculator.release();

%% Plotting BER vs SNR for all channels together

figure(19);
hold on;
semilogy(SNR_dB, BER_AWGN, 'DisplayName', 'AWGN');
semilogy(SNR_dB, BER_ray, 'DisplayName', 'Flat Fast Fading Rayleigh');
semilogy(SNR_dB, BER_ray2, 'DisplayName', 'Frequency-Selective Fast Fading Rayleigh');
semilogy(SNR_dB, BER_ray3, 'DisplayName', 'Flat Slow Fading Rayleigh');
semilogy(SNR_dB, BER_ray4, 'DisplayName', 'Frequency-Selective Slow Fading Rayleigh');
semilogy(SNR_dB, BER_ric, 'DisplayName', 'Flat Fast Fading Rician');
semilogy(SNR_dB, BER_ric2, 'DisplayName', 'Frequency-Selective Fast Fading Rician');
semilogy(SNR_dB, BER_ric3, 'DisplayName', 'Flat Slow Fading Rician');
semilogy(SNR_dB, BER_ric4, 'DisplayName', 'Frequency-Selective Slow Fading Rician');
hold off;
legend('show');
title('BER Comparison For All Channels');
xlabel('SNR in dB');
ylabel('BER');
xlim([1 40]);
grid on;

%% Constellation Diagram of Clean Signal
if exist('modulatedSignal', 'var')
    figure;
    scatterplot(modulatedSignal);
    title('Constellation Diagram of Clean Signal');
else
    disp('Error: modulatedSignal variable is missing for constellation diagram.');
end

%% Constellation Diagram of Received Signal in AWGN Channel
% Replace noisySignal_AWGN with the actual variable from your script
 figure;
    scatterplot(noisySignal);
    title('Constellation Diagram of Received Signal in AWGN Channel');

