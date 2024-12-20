%% Demo of discrete-time systems analysis
%% For the ECE-NTUA Signals&Systems course.

clear all;
close all;

%% Definition of poles/zeros
%% Zeros: 1 (x2), -1
%% Poles: 0, \pm 0.7i

z = [1 1 -1];
p = [0 0.7*i -0.7*i];

z = z(:);
p = p(:);

pause();

%% Plotting the zero-pole diagram

figure
zplane(z,p)

pause();

%% Using zp2tf() to get a [b,a] vector representation
%% so that sum(a_i y[n-i]) = sum(b_i x[n-i)
%% + plotting magnitude/phase response

[b,a] = zp2tf(z,p,1) % arbitrary gain value
figure 
freqz(b,a)

pause();
%% Using the a/b vector representation to compute the impulse and step responses.

figure
impz(b,a)

pause();

figure
stepz(b,a)

pause();

%% Use tf2zp() to revert to the zero-pole representation + cross-check

[z,p,k] = tf2zp(b,a)

pause();
%% Performing filtering (calculate the response to a random input) using 2 methods:
%% a) impulse response calc + convolution
%% b) usage of the built-in filter() function.

t = linspace(0,1,300);
x = 1 + 0.5*sin(10*pi*t) + 0.2*sin(80*pi*t);

figure;
plot(x)

pause();

%a)
h = impz(b,a);
y1 = conv(x,h);

%b)
y2 = filter(b,a,x);

%% Note: filter() outputs a signal of equal length to x
%% impz() + conv() a signal of length = length(x) + length(h) -1

pause();

figure;
plot(y1)
hold on;
plot(y2)

%% How do you explain the filtering result? 