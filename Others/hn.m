%Compute and plot h[n]

W1 = 450;
W2 = 875;
Ws  = 2600;

T = 2*pi/Ws;

N = 79;

range1 = (N-1)/2;
n1 = -range1 : 1 : -1;
h1 = ((1/pi)./n1).*(sin(W1*T.*n1) - sin(W2*T.*n1));
h0 = 1 + 2*(W1 - W2)/Ws;
n2 = 1 : 1 : range1;
h2 = ((1/pi)./n2).*(sin(W1*T.*n2) - sin(W2*T.*n2));
h = [h1,h0,h2];
n = [n1,0,n2];
figure;stem(n,h);grid on;
xlabel('n');ylabel('h[n]');title('Impulse Response of Ideal Bandstop Filter');
