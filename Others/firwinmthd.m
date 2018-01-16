function firwinmthd

% Author - Chamira Edussooriya
% Date - Jan 11, 2017
% Last modified - Apr 18, 2017

clc;
close all;

wc1 = 0.3*pi;
wc2 = 0.5*pi;
wc3 = 0.7*pi;
N = 1000;            % order of the filter
n = 1:N/2;

h = zeros(1,N+1);

h(N/2+1) = (wc1+wc3-wc2)/pi;    % impuse reposne of the multi-band filter
h(N/2+2:end) = (sin(wc1*n)+sin(wc3*n)-sin(wc2*n)) ./ (n*pi);
h(1:N/2) = fliplr(h(N/2+2:end));

wn = zeros(1,N+1);
wn(N/2+1) = 1;

% Windoe function (comment the windows that are not used in the design)

%wn = ones(1,N+1);                           % rectangular window
%wn(N/2+2:end) = 0.5 + 0.5*cos(2*pi*n/N);    % von Hann window
%wn(N/2+2:end) = 0.54 + 0.46*cos(2*pi*n/N);  % Hamming window
%wn(N/2+2:end) = 0.42 + 0.5*cos(2*pi*n/N) + 0.08*cos(4*pi*n/N);  % Blackmann
    % window
wn = chebwin(N+1,200).';                    % Dolph-Chebyshev window
wn(1:N/2) = fliplr(wn(N/2+2:end));

hf = h .* wn;

h
wn
hf

fvtool(hf,1);



