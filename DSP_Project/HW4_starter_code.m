% Code to plot a "continous" signal, sample, downsample, upsample and 
% interpolate (with ideal lowpass filter)
% Taisa Kushner
% Signals and Systems II
% plot a signal and spectrum for HW 4

%% HOW TO RUN A SECTION OF CODE %%
% Sections are separated with '%%' (such as this section is here)
% To "run" a section at a time, click your mouse in the section
% you want to run. This should highlight the code within
% the section, while the code outside of the section remains 
% not highlighted (it is yellow on my computer, but may be different
% depening on your Matlab settings). Then click the "Run Section" button at the top of Matlab
%% 

%% Section of "Base Information" Run this first before any other section 

clear all
close all

% CHANGE THIS PARAMETER fs %
fs = 9000; %sampling frequency (Change this to match what you found)
ts = 1/fs; %sampling interval
% CHANGE THIS PARAMETER fo %
fo = 300; %fundamental freq in hz (Change this to match what you found)
N = fs/fo; %number of samples
to = 1/fo; %fundamental interval 
n = 0:N-1; % the labels of our sequence
nts = n*ts; % the labels of our sampled points in time

M = 4096; %this can be anything, but 4096 is the usual number of adc spots
t = linspace(0,to,M); %sample locations


% CHANGE THIS FUNCTION %
x = cos(600*pi*t)+sin(4000*pi*t)+sin(1800*pi*t); %the "continuous" plot (Change this to match the signal in your homework)

figure(1)
plot(t,x) %plot the "continuous" signal
xlabel('Time (s)') 
ylabel('x(t)')

% CHANGE THIS FUNCTION %
%Here the times n*ts are your sampled location points, eg your "time"
%vector
xs = cos(600*pi*nts)+sin(4000*pi*nts)+sin(1800*pi*nts); %the "continuous" plot (Change this to match the signal in your homework)

hold on
plot(n*ts, xs, 'o') %plot dots at your sampled locations
hold off

display(strcat('The number of samples is : ', string(N)))
%% Run this section to take the discrete time fourier transform of your sampled data
X = fft(xs); %this function takes the fast fourier transform of your sampled data
X = abs(X)/N; %need to divide by the number of samples
X = [fliplr(X(2:end)) X]; % replicate the negative frequencies 
F = -1:1/N:1; %this is just giving us location to plot at

figure(2)
stem(F(2:end-1),X) %stem plots vertical lines (rather than connecting points across the x-axis like plot() does)
title('X|(F)|')
xlabel('F (sample interval in multiples of T)')

%% Run this section to change the labels on your Fourier Transform plot from the prior section
% change the labels on the x-axis from multiples of the time interval to time in seconds by diving by the time interval 
figure(4)
stem(F(2:end-1)/ts,X) %we divide by the time interval so the x-axis is in frequency (hz)
title('X|(F)|')
xlabel('F (sample interval in freq)')


%% Run this section to get the fourier transform of the "continuous" signal
% Note: this is not truly "continuous" but just very high rate of sampling
% a discrete signal, with some cheap tricks to "look" like the continuous spectrum.
% Due to this, you may see some "noise" in the spectrum 
% surrounding your 6 expected frequencies (3 positive, 3 negative). 
% This is due to aliasing of the higher harmonics, you may safely ignore any frequency with 
% power < 0.1, it is noise. 

Tmax = 10*to; %high rate to sample at
tt = linspace(0,Tmax,M);
%CHANGE THIS FUNCTION
xtt = cos(600*pi*tt)+sin(4000*pi*tt)+sin(1800*pi*tt); %Change this function to match the one in your homework, keep the time vector as tt
X = fftshift(fft(xtt)); %fftshift shifts the zero freq component to the center of the plot
X = abs(X)/M;
f = (-M/2:(M/2)-1)/Tmax;
figure(3)
stem(f,X)
xlim([-fs, fs])
title('X|(F)|, power spectrum')
xlabel('frequency, f (Hz)')

%%
%% Run this section to downsample by a factor of SS
%set up vector nans
SS = 1; %Change this to be your downsampling rate (whole number)
Nd = N/SS; % this is our new number of data points
xd = nan(1,N); %create an empty arry of non a number (nan) terms
xd(1:SS:end) = xs(1:SS:end); %pick every SSth value from your original data, and keep it (ignmore the others)

figure(7)
stem(n, xd) %now we plot the data we have kept
xlabel('n')

%drop the terms that are "zero" (in our case it is nan, as it is easier to keep
%track of)
xd(isnan(xd))=[];
%check the length of xd, is it as expected?
display(strcat('length of xd is: ', string(length(xd))))
Nd = length(xd);
%% Run this section ot Upsample by a factor of UU
UU = 3; % Change this value to be your upsampling rate
Nup = Nd*UU;
xup = zeros(1,Nup);
xup(1:UU:end) = xd;
figure(11)
stem(xup)
xlabel('n');
title('up sampled discrete time signal')

%% Run this section to interpolate the down sampled then up sampled data
% this will "fill in" the points that were zero in frequency with values
% based on the sinc function of the surrounding points (in time)

fds = fs/SS*UU; % we get our absolute sampling rate by dividing original sampling rate (Hz) by SS (the value you downsampled at) then mult by UU (rate upsampled at)
tds = 1/fds; % this is the time interval of samples 
fdc = fds/2; %lowpass filter cutoff frequency (1/2 of our sampling rate)
xi = zeros(1,Nup);
ti = 0:tds:(Nup-1)*tds;
arg = fds*ti; %these our our sequence numbers 0 to Nup

% now we "convolve" our sampled points in frequency with the ideal lowpass filter.
% this is the equivalent of multiplying our data points in time with the
% inverse fourier transform of the ideal lowpass filter (the sinc() function)
for m = 0:Nup-1
    xi = xi+xup(m+1)*sinc((arg-2*m*fdc*tds)./UU); %divide by locations where the signal is zero
end

figure(14)
stem(arg*tds, xi)
hold on
%plot the "continuous" signal on top
plot(t,x,'-r')
% ** NOTE ** as the "continuous" signal is not actually continuous
% (discrete points in time) if you are very close to Nyquist sampling rate,
% you may see slight numerical errors between the interpolated signal and
% the "continuous" signal. 

%% Run this section to interopolate the original data
% this is a bit lame, as we just get the same data points that we sampled at (not filling
% anything in)
fc = fs/2; %/2; %lowpass filter cutoff frequency
xi = zeros(1,N);
ti = 0:ts:(N-1)*ts;
arg = fs*ti;
for m = 0:N-1
    xi = xi+xs(m+1)*sinc(arg-2*m*fc*ts);
end
figure(12)
stem(arg*ts, xi)
hold on

%note how the sinc(arg-2*n*fc*ts) is 1 at all times

% plot the "continous" function too
plot(t,x,'-r')

%% Run this section to try and convolve with a lowpass filter of the wrong cutoff frequency
% This is for fun. You can see how the signal we get back in time now has frequencies present that should
% not be. Notice how the shape changes

fc = 6000; % let all the frequencies up to our sampling rate through! (make this anything)
xi = zeros(1,N);
ti = 0:ts:(N-1)*ts;
arg = fs*ti;

for m = 0:N-1
    xi = xi+xs(m+1)*sinc((arg-2*m*fc*ts)); %sinc function centered at the m-th sample point
end

figure(12)
stem(arg*ts, xi, '-o')
hold on
%plot the continous function too
plot(t,x,'-r')


%%