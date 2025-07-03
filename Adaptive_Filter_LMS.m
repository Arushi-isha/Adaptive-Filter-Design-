clc;                % Clear the command window
clear all;          % Remove all variables from workspace
close all;          % Close all figure windows

% ========================
% Signal and Sampling Setup
% ========================
fs = 500;           % Sampling frequency in Hz
f = 5;              % Frequency of clean sine wave in Hz
n = 1000;           % Total number of samples
t = (0:n-1) / fs;   % Time vector corresponding to sampling rate

% ========================
% Generate Clean and Noisy Signals
% ========================
clean_signal = sin(2 * pi * f * t);       % Generate a clean sine wave signal
noise = 0.5 * randn(1, n);                % Generate white Gaussian noise with std dev = 0.5
noisy_signal = clean_signal + noise;     % Add noise to clean signal to get the observed (noisy) signal

desired_signal = clean_signal;            % The desired signal is the original clean signal

% =========================
% LMS Adaptive Filter Initialization
% =========================
mu = input("Enter the step size: ");         % Step size (learning rate) for LMS algorithm
M = 64;                                      % Filter order (number of coefficients/taps)
w = zeros(1, M);                             % Initialize filter weights to zero
output_lms = zeros(1, n);                    % Initialize output of LMS filter
error_lms = zeros(1, n);                     % Initialize error signal
SqErr = zeros(1, n);                         % Squared error at each time step

% =========================
% LMS Adaptive Filtering Process
% =========================
for i = M+1:n
    x = noisy_signal(i-M:i-1);                     % Extract a segment (window) of the noisy input signal (length M)
    y = w * x.';                                   % Calculate filter output using current weights
    error_lms(i) = desired_signal(i) - y;          % Calculate error between desired and output
    w = w + 2 * mu * error_lms(i) * x;             % Update weights using LMS update rule
    output_lms(i) = y;                             % Store the filter output
end

% =========================
% FIR Filter Design using Window Methods
% =========================
cutoff = 10;      % Cutoff frequency for FIR filter (Hz)
numtaps = 64;     % Number of filter coefficients (taps)

% FIR filter using Hamming window
fir_hamming = fir1(numtaps-1, cutoff/(fs/2), 'low', hamming(numtaps));
output_hamming = filter(fir_hamming, 1, noisy_signal);   % Apply filter to noisy signal

% FIR filter using Blackman window
fir_blackman = fir1(numtaps-1, cutoff/(fs/2), 'low', blackman(numtaps));
output_blackman = filter(fir_blackman, 1, noisy_signal); % Apply filter to noisy signal

% FIR filter using Rectangular window
fir_rect = fir1(numtaps-1, cutoff/(fs/2), 'low', rectwin(numtaps));
output_rect = filter(fir_rect, 1, noisy_signal);         % Apply filter to noisy signal

% =========================
% Plot Comparison: LMS vs Hamming Window FIR
% =========================
figure('Name','LMS vs FIR (Hamming)','NumberTitle','off');
subplot(3,1,1);
plot(t, clean_signal, 'g', 'LineWidth', 1.5);          % Plot clean signal
hold on;
plot(t, noisy_signal, 'r', 'LineWidth', 0.7);          % Plot noisy signal
title('Clean Signal vs Noisy Signal');                 % Title and legend
legend('Clean Signal', 'Noisy Signal','Location','best');
grid on;

subplot(3,1,2);
plot(t, clean_signal, 'g', t, output_lms, 'b','LineWidth', 1.5);  % Plot LMS output
title('LMS Adaptive Filter Output');
legend('Clean Signal', 'LMS Output','Location','best');
grid on;

subplot(3,1,3);
plot(t, clean_signal, 'g', t, output_hamming, 'm','LineWidth', 1.5);  % Plot Hamming FIR output
title('Static FIR (Hamming Window)');
legend('Clean Signal', 'Hamming Output','Location','best');
grid on;

% =========================
% Plot Comparison: LMS vs Blackman Window FIR
% =========================
figure('Name','LMS vs FIR (Blackman)','NumberTitle','off');
subplot(3,1,1);
plot(t, clean_signal, 'g', 'LineWidth', 1.5);
hold on;
plot(t, noisy_signal, 'r', 'LineWidth', 0.7);
title('Clean Signal vs Noisy Signal');
legend('Clean Signal', 'Noisy Signal','Location','best');
grid on;

subplot(3,1,2);
plot(t, clean_signal, 'g', t, output_lms, 'b','LineWidth', 1.5);
title('LMS Adaptive Filter Output');
legend('Clean Signal', 'LMS Output','Location','best');
grid on;

subplot(3,1,3);
plot(t, clean_signal, 'g', t, output_blackman, 'c','LineWidth', 1.5);
title('Static FIR (Blackman Window)');
legend('Clean Signal', 'Blackman Output','Location','best');
grid on;

% =========================
% Plot Comparison: LMS vs Rectangular Window FIR
% =========================
figure('Name','LMS vs FIR (Rectangular)','NumberTitle','off');
subplot(3,1,1);
plot(t, clean_signal, 'g', 'LineWidth', 1.5);
hold on;
plot(t, noisy_signal, 'r', 'LineWidth', 0.7);
title('Clean Signal vs Noisy Signal');
legend('Clean Signal', 'Noisy Signal','Location','best');
grid on;

subplot(3,1,2);
plot(t, clean_signal, 'g', t, output_lms, 'b','LineWidth', 1.5);
title('LMS Adaptive Filter Output');
legend('Clean Signal', 'LMS Output','Location','best');
grid on;

subplot(3,1,3);
plot(t, clean_signal, 'g', t, output_rect, 'k','LineWidth', 1.5);
title('Static FIR (Rectangular Window)');
legend('Clean Signal', 'Rectangular Output','Location','best');
grid on;

% =========================
% Convergence Curve (MSE in dB)
% =========================
eps_val = 1e-12;  % Small value to avoid log(0)

% Compute cumulative moving average of squared error to show convergence over time
SqErr = error_lms.^2; 
MSE_cumulative = cumsum(SqErr) ./ (1:n);
MSE_dB = 10 * log10(MSE_cumulative + eps_val);

% Plot convergence curve
figure('Name','LMS Convergence Curve','NumberTitle','off');
plot(MSE_dB, 'LineWidth', 2);
xlabel('Number of Iterations');
ylabel('MSE (dB)');
title('Convergence of LMS Adaptive Filter (MSE in dB)');
grid on;
