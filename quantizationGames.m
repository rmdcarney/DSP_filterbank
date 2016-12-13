%% Full filter bank


%%% Toy Data %%%
 x = linspace(0,1,2048); %Create 1000 samples from [0,pi/2)
 %y = 0.5*(sin(200*pi*x) + sin(1638*pi*x)); %Make an input dataset with 2 distinct peaks
% freq = linspace(-0.5,0.5,2048); %Only want to plot 0-0.5
[y, Fs] = audioread('orinoccio.wav');
%[y, Fs] = audioread('thank.wav');

b_h0 = 0.125*[-1 2 6 2 -1];
a_h0 = 1;
b_g0 = 0.5*[1 2 1];
a_g0 = 1;  

b_h1 = 0.5*[1 -2 1];
a_h1 = 1;
b_g1 = 0.125*[1 2 -6 2 1];
a_g1 = 1;


%==== ANALYSIS FILTER BANK ========= 

%---- High pass branch ----

%Use matlab's filter-maker to high pass filter
U1 = filter(b_h1, a_h1, y);

%Downsample high pass branch
Y1 = U1(1:2:end);

%---- Low pass branch ----

%Use matlab's filter-maker to low pass filter
U0 = filter(b_h0, a_h0, y);

%Downsample low pass branch
Q0 = U0(1:2:end);

%Filter again
W0 = filter(b_h0,a_h0,Q0);
W1 = filter(b_h1,a_h0,Q0);

%Downsample each subbranch
Y00 = W0(1:2:end);
Y01 = W1(1:2:end);

%==== Quantization =====

%--- High band ----
bits_Y1 = 1;
scale_Y1 = max(abs(Y1))/(1-pow2(-bits_Y1));
Y1_Q = scale_Y1*double(fixed(bits_Y1,Y1/scale_Y1));

%--- Low-high band ---
bits_Y01 = 6;
scale_Y01 = max(abs(Y01))/(1-pow2(-bits_Y01));
Y01_Q = scale_Y01*double(fixed(bits_Y01,Y01/scale_Y01));

%--- Low-low band ---
bits_Y00 =8;
scale_Y00 = max(abs(Y00))/(1-pow2(-bits_Y00));
Y00_Q = scale_Y00*double(fixed(bits_Y00,Y00/scale_Y00));
% 
% Y1_Q = Y1;
% Y00_Q = Y00;
% Y01_Q = Y01;

%======= SYNTHESIS FILTER BANK =======

%---- High pass branch ----

%Delay the high pass branch
T1 = filter([0 0 0 1],1,Y1_Q');
T1 = T1';
%Upsample the high pass branch
R1 = zeros(1,length(T1)*2);
R1(1:2:end) = T1(1:1:end);

%Filter the upsampled signal in high pass branch
S1 = filter(b_g1,1,R1);

%---- Low pass branch ----

%Upsample incoming signals
P01 = zeros(1,length(Y01_Q)*2);
P01(1:2:end) = Y01_Q(1:1:end);

P00 = zeros(1,length(Y00_Q)*2);
P00(1:2:end) = Y00_Q(1:1:end);

%Filter upsampled signal to interpolate
K01 = filter(b_g1,1,P01);
K00 = filter(b_g0,1,P00);

%Sum the filtered signals
K0 = K00 + K01;

%Upsample summed branch
R0 = zeros(1,length(K0)*2);
R0(1:2:end) = K0(1:1:end);

%Filter 
S0 = filter(b_g0,1,R0);

%---- Sum two branches ----
V = S0 + S1;

expNoise = mean((V(10:end)-y(1:end-9)').^2);
expY = mean(y.^2);
SQNR = expY/expNoise
%% Plotting %%

% % Input
% [H, w] = freqz(y);
% 
% %Original
% figure(1);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of orinoccio audio sample');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% %legend('0.5 (sin(200\pix) + sin(1638\pix))');
% saveas(figure(1),'orinoccio_original.eps','epsc');
% 
% % ========= HIGH PASS ANALYSIS ==============
% 
% % High pass filter
% [H, w] = freqz(U1);
% 
% figure(2);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal HP filtered with H_1(z)');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(2),'testSignal_U1.eps','epsc');
% 
% % High pass downsample
% [H, w] = freqz(Y1);
% 
% figure(3);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal HP filtered with H_1(z) and downsampled by 2');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(3),'testSignal_Y1.eps','epsc');
% 
% % ========= LOW PASS ANALYSIS ==============
% 
% % Low pass filter
% [H, w] = freqz(U0);
% 
% figure(4);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal LP filtered with H_0(z)');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(4),'testSignal_U0.eps','epsc');
% 
% % Decimation
% [H, w] = freqz(Q0);
% 
% figure(5);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal LP filtered with H_0(z) and downsampled by 2');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(5),'testSignal_Q0.eps','epsc');
% 
% % High pass Low pass filter
% [H, w] = freqz(W1);
% 
% figure(6);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal HP filtered H_1(z) in second stage');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(6),'testSignal_W1.eps','epsc');
% 
% % High pass Low pass downsampled 
% [H, w] = freqz(Y01);
% 
% figure(7);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal downsampled by 2 in second stage HP branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(7),'testSignal_Y01.eps','epsc');
% 
% % Low pass Low pass filter
% [H, w] = freqz(W0);
% 
% figure(8);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal LP filtered H_0(z) in second stage');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(8),'testSignal_W0.eps','epsc');
% 
% % Low pass Low pass downsampled 
% [H, w] = freqz(Y00);
% 
% figure(9);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal downsampled by 2 in second stage LP branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(9),'testSignal_Y00.eps','epsc');
% 
% % ========= HIGH PASS SYNTHESIS ==============
% 
% %Upsampled
% [H, w] = freqz(R1);
% 
% figure(10);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal upsampled by 2, HP branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(10),'testSignal_R1.eps','epsc');
% 
% % HP filtered
% [H, w] = freqz(S1);
% 
% figure(11);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal HP filtered G_1(z) in sythesis branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(11),'testSignal_S1.eps','epsc');
% 
% % ========= LOW PASS Synthesis ==============
% 
% %Upsampled
% [H, w] = freqz(P01);
% 
% figure(12);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal upsampled by 2, 2nd level synthesis branch, HP');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(12),'testSignal_P01.eps','epsc');
% 
% %HP Filtered
% [H, w] = freqz(K01);
% 
% figure(13);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal HP filtered G_1(z) in sythesis branch, 2nd level');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(13),'testSignal_K01.eps','epsc');
% 
% %Upsampled
% [H, w] = freqz(P00);
% 
% figure(14);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal upsampled by 2, 2nd level synthesis branch, LP');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(14),'testSignal_P00.eps','epsc');
% 
% %HP Filtered
% [H, w] = freqz(K00);
% 
% figure(15);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal LP filtered G_0(z) in sythesis branch, 2nd level');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(15),'testSignal_K00.eps','epsc');
% 
% %Combined
% [H, w] = freqz(K0);
% 
% figure(16);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal in 2nd level summed');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(16),'testSignal_K0.eps','epsc');
% 
% %Upsampled
% [H, w] = freqz(R0);
% 
% figure(17);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal upsampled by 2 - LP synthesis branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(17),'testSignal_R0.eps','epsc');
% 
% %LP Filtered
% [H, w] = freqz(S0);
% 
% figure(18);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of test signal LP filtered G_0(z) in sythesis branch');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(18),'testSignal_S0.eps','epsc');
% 
% %Final Reconstructed Signal
% [H, w] = freqz(V);
% 
% figure(19);
% plot(w/(2*pi),20*log10(abs(H)));
% grid on;
% title('PSD of reconstructed orinoccio file');
% xlabel('Normalized Frequency (\nu) [2\pi rad/sample]')
% ylabel('Magnitude [dB]');
% saveas(figure(19),'thank_V.eps','epsc');

