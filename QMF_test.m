%% Mini filterbank test

%%% Toy Data %%%
x = linspace(0,1,2048); %Create 1000 samples from [0,pi/2)
y = 0.5*(sin(200*pi*x));% + sin(1638*pi*x)); %Make an input dataset with 2 distinct peaks

%==== ANALYSIS FILTER BANK ========= 

%---- Low pass branch ----

%Use matlab's filter-maker to low pass filter
b_h0 = 0.125*[-1 2 6 2 -1];
a_h0 = 1;
A0 = filter(b_h0,a_h0,y);

%Downsample
B0 = A0(1:2:end);

%Upsample
C0 = zeros(1,length(B0)*2);
C0(1:2:end) = B0(1:1:end);

%Low pass filter again
b_g0 = 0.5*[1 2 1];
a_g0 = 1;
D0 = filter(b_g0,a_g0,C0);

%---- High pass branch ----

%High pass filter signal
b_h1 = 0.5*[1 -2 1];
a_h1 = 1;
A1 = filter(b_h1,a_h1,y);

%Downsample
B1 = A1(1:2:end);

%Upsample
C1 = zeros(1,length(B1)*2);
C1(1:2:end) = B1(1:1:end);

%High pass filter again
b_g1 = 0.125*[1 2 -6 2 1];
a_g1 = 1;
D1 = filter(b_g1,a_g1,C1);

%---- Sum ----
Z = D0 + D1;