function BandStopAba

%Author: Abarajithan G. : 150001C

clc;
close all;
clear all;

%Parameters
Ap_ = 0.05;     % Min Passband Ripple (db)
Aa_ = 40;       % Min Stopband Attenuation (db)
Wp1 = 400;      % Lower passband Freq (rad/s)
Wa1 = 500;      % Lower stopband Freq (rad/s)
Wa2 = 800;      % Higher stopband Freq (rad/s)
Wp2 = 950;      % Higher stopband Freq (rad/s)
Ws  = 2600;     % Sampling Freq (rad/s)

Bt = min(Wa1-Wp1, Wp2-Wa2);

% 2. Choose Delta

delta_a = 10^(-Aa_/20);
c = 10 ^ ( Ap_ / 20);
delta_p = (c-1)/(c+1);
delta = min(delta_a, delta_p);

% 3. Get Aa from delta
Aa = -20*log10(delta);            % Actual stopband attenuation

% 4. Calculate alpha
if (Aa <= 21)
    alpha = 0;
elseif ((21 < Aa) && (Aa <= 50))
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa - 8.7);
end

% 5. Calculate D and N

if (Aa <= 21)
    D = 0.9222;
else
    D = (Aa - 7.95)/14.36;
end

N = ceil ( Ws * D / Bt +1);

if (mod(N,2) == 0)
    N = N + 1;
end

% 6. Form Kaiser window

n  = -(N-1)/2 :1: (N-1)/2 ;

beta = alpha*(1 - (2*n/(N-1)).^2).^0.5;
Ibeta = 0; 
Ialpha = 0;

for k = 1 : 15						%10-15 are enough
    Ibeta = Ibeta + ((1/factorial(k))*(beta/2).^k).^2;
    Ialpha = Ialpha + ((1/factorial(k))*(alpha/2)^k)^2;
end
Ibeta = Ibeta + ones(1,numel(Ibeta));
Ialpha = Ialpha + ones(1,numel(Ialpha));

wk = Ibeta ./ Ialpha;


figure;
stem(n,wk);
xlabel('n');
ylabel('w_k[n]');
title('Kaiser Window');

% 7. Impulse response h[n] of Bandstop

W1 = (Wp1+Wa1)/2;
W2 = (Wa2+Wp2)/2;

n1 = -(N-1)/2 : 1 : -1; %Negative range
n2 = 1 : 1 : (N-1)/2;   %Positive range

h1 = ((1/pi)./n1).*(sin((2*pi)*W1/Ws.*n1) - sin((2*pi)*W2/Ws.*n1));
h0 = 1 + 2*(W1 - W2)/Ws;
h2 = ((1/pi)./n2).*(sin((2*pi)*W1/Ws.*n2) - sin((2*pi)*W2/Ws.*n2));

h_ideal = [h1,h0,h2];
n = [n1,0,n2];

figure;
stem(n,h_ideal);

xlabel('n');
ylabel('h[n]');
title('Impulse Response of Expected Ideal Filter');

% Applying the Window

h_final = h_ideal .* wk;    % Non causal Filter

figure;
stem(n,h_final);
xlabel('n');
ylabel('h[n]');
title('Non-causal Filter: Impulse Response');


stem( 1:size(h_final, 2) ,h_final);

xlabel('n');
ylabel('h[n]');
title('Causal Filter: Impulse Response')

% Inspect the filter's magnitude, phase response, group delay and phase
% delay

fvtool(h_final);

% Creating Given Excitation

n = 0 : 1 : 250;
l = size(n,2);

We1 = W1/2
We2 = (W2+W1)/2
We3 = (Ws/2+W2)/2

T = 2*pi/Ws;
x = sin(We1*T*n) + sin(We2*T*n) + sin(We3*T*n); %Excitation Signal

% Taking Discrete Fourier Transform of Excitation

NFFT = 2^nextpow2(l); 
X = fft( x , NFFT) / l;
f = (Ws)/2 * linspace(0,  1, NFFT/2+1);
mag = 2*abs(X(1:NFFT/2+1));

figure;
plot(f , mag);

title('DFT of Excitation Signal')
xlabel('Frequency (rad/s)')
ylabel('|X(w)|');

%   Applying the filter to the excitation

x_filtered = conv(x , h_final);
X_filtered = fft(x_filtered,NFFT)/l;
mag = 2*abs(X_filtered(1:NFFT/2+1));

figure;
plot(f,mag);
title('DFT of Excitation Passed Through Filter')
xlabel('Frequency (rad/s)')
ylabel('|X_k(w)|');

%Passing the excitation through an ideal bandstop filter

x_idealFiltered = sin(We1*T*n) + sin(We3*T*n);
X_idealFiltered = fft(x_idealFiltered,NFFT)/l;
mag = 2*abs(X_idealFiltered(1:NFFT/2+1));

figure;
plot(f,mag);
title('DFT of Excitation through an Ideal Filter')
xlabel('Frequency (rad/s)')
ylabel('|X_i(w)|');


%   Plotting Excitation

figure;
plot(x);
title('Time doman Response of Excitation')
xlabel('Time(s)')
ylabel('x[n]');
axis([50,250, -2.5,2.5]);

figure;
plot(x_filtered);
title('Time doman Response of Filtered Excitation')
xlabel('Time(s)')
ylabel('x[n]');
axis([50,250, -2.5,2.5]);


alpha
D
N
Bt












    