%% Full filter bank


%%% Toy Data %%%
x = linspace(0,1,2048); %Create 1000 samples from [0,pi/2)
y = 0.5*(sin(200*pi*x));% + sin(1638*pi*x)); %Make an input dataset with 2 distinct peaks
freq = linspace(-0.5,0.5,2048); %Only want to plot 0-0.5

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

%======= SYNTHESIS FILTER BANK =======

%---- High pass branch ----

%Delay the high pass branch
T1 = [[0 0 0],Y1(1:end-3)];

%Upsample the high pass branch
R1 = zeros(1,length(T1)*2);
R1(1:2:end) = T1(1:1:end);

%Filter the upsampled signal in high pass branch
S1 = filter(b_g1,1,R1);

%---- Low pass branch ----

%Upsample incoming signals
P01 = zeros(1,length(Y01)*2);
P01(1:2:end) = Y01(1:1:end);

P00 = zeros(1,length(Y00)*2);
P00(1:2:end) = Y00(1:1:end);

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


