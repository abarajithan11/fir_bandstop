
t = 0:1/50:10-1/50;                     
x = sin(2*pi*15*t) + sin(2*pi*20*t);
figure;
plot(t,x)

y = fft(x);     
f = (0: length(y)-1)  *50 / length(y);

figure;
plot(f,abs(y))
title('Magnitude')

n = length(x);                         
fshift = (-n/2:n/2-1)*(50/n);
yshift = fftshift(y);
figure;
plot(fshift,abs(yshift));