% Anthony Bruno
% ECE 316 Plot
% HW #1 - Problem 3 plot of x(t)

clear
clc

t = -2*pi:pi/100:2*pi;
x = (8/(2*pi))*(sin(2*pi.*t).*heaviside(t) - sin(2*pi.*(t-1)).*heaviside(t-1));

figure
plot(t,x) 
grid on
title('Problem 3 - Plot of x(t)')
xlabel('t')
ylabel('x(t)')