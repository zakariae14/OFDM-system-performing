clc;
clear;
close all;
% Initializing
data_bits = 128;
nsc = 4;
sub_car = 4;
size_block = 32;
length_cp = floor(0.15 * size_block);
N = 1000;  % Number of iteration %% to have more smoothly, increase N


% AWGN Channel
% Channel parameters
SNRdB = -20:1:5 ; % signal-to-noise ratio in dB
SER = [];
SER_per_SNR = 0;
BER = [];
BER_per_SNR = 0;

% Transmitter %

% Generate random data source to be transmitted of size (128 , N)

bits = randsrc(1,log2(nsc)*N*data_bits, 0:1);
data_src = zeros(1,N*data_bits);
up = 1;
for b=1:log2(nsc)*N*data_bits
    if mod(b,2) == 0
        symbol = 2*bits(b-1) + bits(b);
        data_src(up) = symbol;
        up = up + 1;
    end  
end


for f=1:length(SNRdB)
    qpsk_demodulated_data_list = [];


for e=1:N
qpskmod_data = pskmod(data_src(data_bits*(e-1) +1:data_bits*e), nsc,pi/nsc);

% Conversion of series data stream into four parallel data stream to form four sub carriers
series_to_parallel = reshape(qpskmod_data, data_bits/sub_car,sub_car);

% Inverse Fast Fourier Transform of Sub-Carriers
start_cp = size_block-length_cp;

% Adding Cyclic Prefix
sc_ifft = zeros (size_block, sub_car);
cyclic_prefix = zeros (length_cp, sub_car);
Append_prefix = zeros (size_block + length_cp, sub_car);
for i=1:sub_car
%sc_ifft(:,i) = ifft((series_to_parallel(:,i)),size_block); % ifft
sc_ifft(:,i) = ifrft((series_to_parallel(:,i)),size_block); % ifrft = inverse fast fractional Fourier transform
for j=1:length_cp
cyclic_prefix(j,i) = sc_ifft(j+start_cp,i);
end
Append_prefix(:,i) = vertcat( cyclic_prefix(:,i), sc_ifft(:,i));
end


%Conversion to serial for transmission
[rows_Append_prefix,cols_Append_prefix]=size(Append_prefix);
ofdm_datalen = rows_Append_prefix*cols_Append_prefix;

% Trasmitting of OFDM Signal
ofdm_signal = reshape(Append_prefix, 1, ofdm_datalen);


% Calculate noise power spectral density
SNR = 10.^(SNRdB(f)/10); % convert SNR from dB to linear scale

sigPower = mean(abs(ofdm_signal).^2); % calculate signal power
noisePower = sigPower / SNR; % calculate noise power
N0 = noisePower / 2; % calculate one-sided noise power spectral density

% Generate complex noise with N0 power spectral density
noise = sqrt(N0/2)*(randn(size(ofdm_signal)) + 1i*randn(size(ofdm_signal)));

% Add noise to OFDM signal
received_sig = ofdm_signal + noise;

           % Reciever %
received_sig_parallel = reshape(received_sig,rows_Append_prefix, cols_Append_prefix);

% Removing Cyclic Prefix
received_sig_parallel(1:length_cp,:)=[];


% Recieved signal FFT / FRTF
data_fft = zeros (size_block, 1);
for j=1:sub_car
%data_fft(:,j) = fft(received_sig_parallel(:,j),size_block); % fft
data_fft(:,j) = frft(received_sig_parallel(:,j),size_block); % frft = fast fractional Fourier transform
end


% Reconstruced Signal
% Converting to Serial and Demodulation
received_ser_data = reshape(data_fft, 1,[]);
qpsk_demodulated_data = pskdemod(received_ser_data, nsc,pi/nsc);


qpsk_demodulated_data_list = cat(2,qpsk_demodulated_data_list, qpsk_demodulated_data);
end

bits_demod = reshape(cat(1,fix(qpsk_demodulated_data_list/2),mod(qpsk_demodulated_data_list,2)),1,2*length(qpsk_demodulated_data_list));
%for s=1:length(data_src)
 %   symbole = qpsk_demodulated_data_list(s);
  %  bits_demod = cat(2,bits_demod,[fix(symbole/2),mod(symbole,2)]);
%end

% Calculate SER
num_errors = sum(xor(data_src, qpsk_demodulated_data_list));
SER_per_SNR = num_errors/length(data_src);
SER = cat(2,SER,SER_per_SNR);

% Calculate BER
num_errors = sum(abs(bits-bits_demod));
BER_per_SNR = num_errors/length(bits);
BER = cat(2,BER,BER_per_SNR);
end


figure;
semilogy(SNRdB, SER,'-o','Linewidth',1);
hold on;
semilogy(SNRdB, BER,'-o','Linewidth',1);
xlabel('SNR (dB)');
ylabel('Symbol/Bit Error Rate');
legend('SER','BER')
title('SER vs BER for QPSK in OFDM (4 subcarriers)');
grid;
