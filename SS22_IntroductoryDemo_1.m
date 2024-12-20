%% Signal Processing using MATLAB, Part 1
% Brief demonstration of various signal processing techniques using MATLAB
% Demo Subjects: 
%    1. Introduction to MATLAB
%    2. Introduction to "Signal and Systems" basic MATLAB functions 
%     
% Reference: Signals & Systems, (Openheim and Whilsky)
%            Signal Processing First (McCellan, Schafer, Yoder) 
% Project: du
% Course: Signals & Systems, ECE, NTUA
% http://cvsp.cs.ntua.gr/courses/systems
%
% Created: July 01, 2007, Updated: May 15, 2009, Dimitriadis Dimitris (ddim@cs.ntua.gr)Dimitriadis Dimitris
% Updated: 22 April 2010, Georgios Evangelopoulos (gevag@cs.ntua.gr)
% Updated: 11 Nov., 2018, Nancy Zlatintsi (nzlat@cs.ntua.gr)
% Updated: 14 Dec., 2020, Nancy Zlatintsi (nzlat@cs.ntua.gr)


%
%run
%clear all; close all; clc;
%cd ('C:\Users\nancy\Documents\work\courses\ss-2018\gevag\ss2010\demo matlab\demo');

%run('SS20_IntroductoryDemo_1')
% SS20_IntroductoryDemo_1
disp('Welcome to MATLAB ')
disp('For any help, TYPE "help" ')
disp(' ')
pause

% ctrl + c 
%!
disp('plot a figure and listen to the sound')
disp('[yc,fs] = audioread(''piano_E4.wav'')')

[yc,fs] = audioread('piano_E4.wav');
figure
plot(yc);
sound(yc,fs);
axis('tight')
pause
disp(' ')
disp(' ')

disp('plot another figure and listen to the sound')
disp(' [yc,fs] = audioread(''viola_note.wav'') ')
[yc,fs] = audioread('viola_note.wav');
figure
plot(yc);
sound(yc,fs);

pause
disp(' ')
disp(' ')

% i.e., help command
disp('how to use help command and make your scripts reusable by others')
disp(' ')
disp('help(''SS20_IntroductoryDemo_1'')')
help('SS20_IntroductoryDemo_1')
disp(' ')
pause
disp(' ')
disp(' ')

disp('help(''conv'')')
help('conv')
disp(' ')
pause
disp(' ')
disp(' ')

disp('pi = ')
pi
disp(' ')


disp('3')
3
disp(' ')
pause
disp(' ')
disp(' ')

disp('eps')
eeps = eps
disp(' ')
pause
disp(' ')
disp(' ')


disp('iinf = 1/0')
iinf = 1/0
disp(' ')
pause
disp(' ')
disp(' ')

disp('naan = 0/0')
naan = 0/0
disp(' ')
pause
disp(' ')
disp(' ')
%% Introduction to Vectors, Matrices
% This slide serves as an introduction to the basic MATLAB
% definitions of matrices and vectors.
% 
% We show how to define a vector and it's transpose. Then, we define a
% matrix setting it's elements or by concatenating different vectors.
% 
% We introduce matrix definitions like: ones, zeros, and all the available operations 


disp('Welcome to MATLAB ')
disp('For any help, TYPE "help" ')
disp(' ')
pause

% Definition of Vectors and Matrices
disp('Vector = [1 3 4 4 12 22 15 33 3 3 4 12]')
Vector = [1 3 4 4 12 22 15 33 3 3 4 12]
pause
disp(' ')
disp(' ')

% Definition of Vectors and Matrices
disp('alternatively using commas')
disp('Vector = [1, 3, 4, 4, 12, 22, 15, 33, 3, 3, 4, 12]')
Vector = [1, 3, 4, 4, 12, 22, 15, 33, 3, 3, 4, 12]
pause
disp(' ')
disp(' ')

disp('transpose vector using ''''')
disp('TransVector = [1 3 4 4 12 22 15 33 3 3 4 12]'' ')
TransVector = [1 3 4 4 12 22 15 33 3 3 4 12]'
pause
disp(' ')
disp(' ')

disp('transpose vector using ; ')
disp('TransVector = [1; 3; 4; 4; 12; 22; 15; 33; 3; 3; 4; 12] ')
TransVector = [1; 3; 4; 4; 12; 22; 15; 33; 3; 3; 4; 12]
pause
disp(' ')
disp(' ')


disp('MATRIX creation')
disp(' ')
disp(' ')
disp('Mat1a = [1 2 3; 33 3 23; 1 11 23]'ra)
Mat1a = [1 2 3; 33 3 23; 1 11 23]
pause
disp(' ')
disp(' ')


disp('Mat1b = [[1 2 3]'', [33 3 23]'', [1 11 23]''] ')
Mat1b = [[1 2 3]', [33 3 23]', [1 11 23]']
pause
disp(' ')
disp(' ')


disp('MATRIX concatenation')
disp('Mat2a = [Mat1a, Mat1a]')
Mat2a = [Mat1a, Mat1a]
pause 
disp(' ')
disp(' ')


disp('Mat2b = [Mat1a; Mat1a]')
Mat2b = [Mat1a; Mat1a]
pause 
disp(' ')
disp(' ')

% predifined commands
disp('Predifined commands')
disp(' ')
disp(' ')

disp('I=ones(4,3)')
I=ones(4,3)
pause
disp(' ')
disp(' ')


disp('Z=zeros(5,1)')
Z=zeros(5,1)
pause
disp(' ')
disp(' ')


disp('Z=zeros(1,5)')
Z=zeros(1,5)
pause
disp(' ')
disp(' ')


%Vector adressing
disp('Vector Adressing')
disp(' ')
disp(' ')

disp('Vector = [1 3 4 4 12 22 15 33 3 3 4 12]')
Vector = [1 3 4 4 12 22 15 33 3 3 4 12]
pause
disp(' ')
disp(' ')

disp('x1 = Vector(4)')
x1 = Vector(4)
pause
disp(' ')
disp(' ')

disp('x2 = Vector(1:2:end)')
x2 = Vector(1:2:end)
pause
disp(' ')
disp(' ')


disp('len_sig = length(x2)')
len_sig = length(x2)
pause
disp(' ')
disp(' ')

disp('lin_ex = linspace(1,10,5)')
lin_ex = linspace(1,10,5)
pause
disp(' ')
disp(' ')

disp('lin_ex = linspace(1,10,10)')
lin_ex = linspace(1,10,10)
pause
disp(' ')
disp(' ')


%Elementary Vector and Matrix Operations
disp('Elementary Vector and Matrix Operations')
disp(' ')
disp('I+I')
I+I
pause
disp(' ')
disp(' ')


disp('I-I')
I-I
pause
disp(' ')
disp(' ')


disp('Mat2a*Mat2b')
Mat2a*Mat2b
pause
disp(' ')
disp(' ')


disp('Mat2b*Mat2a')
Mat2b*Mat2a
pause
disp(' ')
disp(' ')


% Pointwise Operations
disp('Mat2a.^2')
Mat2a.^2
pause
disp(' ')
disp(' ')

disp('Mat2a.*Mat2a')
Mat2a.*Mat2a
pause
disp(' ')
disp(' ')

disp('Vector.^2= [1 3 4 4 12 22 15 33 3 3 4 12].^2')
Vector = [1 3 4 4 12 22 15 33 3 3 4 12]
Vector.^2
pause
disp(' ')
disp(' ')


disp('TransVector.^2')
TransVector.^2
pause

disp('sqrt(TransVector.^2)')
sqrt(TransVector.^2)
pause
disp(' ')
disp(' ')


%% Introduction To Basic Function and Control Flows
% Introducing the "plot" tool, for-loops and IF-THEN commands
%
% Plot is a powerfull tool and the begginers are strongly urgued to
% read the manual. To do so, type "help plot" and "doc plot"
%
% Introduction to Matlab control-flow commands and options
%

disp('Plotting!!!')
disp(' ')
disp(' ')
Vector = [1 3 4 4 12 22 15 33 3 3 4 12];
clf;

%Plotting Vectors
disp('plot(Vector); ')
disp(' ')
disp('title(''Plot of a Vector''); ')
disp('xlabel(''Index k''); ')
disp('ylabel(''Amplitude Vector[k]''); ')
disp('shg')
plot(Vector);
title('Plot of a Vector');
xlabel('Index k');
ylabel('Amplitude of Vector[k]');
shg % makes the current figure visible and places it in front of all other figures on the screen.
pause
disp(' ')
disp(' ')


disp('plot( Vector'',''r''); ')
disp('title(''Plot of the Transposed Vector''); ')
plot( Vector','r');
title('Plot of the Transposed Vector');
shg
pause
disp(' ')
disp(' ')


disp('plot(Vector'',''rd-''); ')
disp('title(''Plot of a Transposed Vector''); ')
plot(Vector','rd-');
title('Plot of a Transposed Vector');
shg
pause
disp(' ')
disp(' ')
disp('Parameters of the Plot function')
disp('            b     blue                 .     point                -     solid   ')
disp('            g     green                o     circle               :     dotted  ')
disp('            r     red                  x     x-mark               -.    dashdot ') 
disp('            c     cyan                 +     plus                 --    dashed  ')   
disp('            m     magenta              *     star               (none)  no line ')
disp('            y     yellow               s     square ')
disp('            k     black                d     diamond ')
disp('                                       v     triangle (down) ')
disp('                                       ^     triangle (up) ')
disp('                                       <     triangle (left) ')
disp('                                       >     triangle (right) ')
disp('                                       p     pentagram ')
disp('                                       h     hexagram ')
disp(' ')
disp(' ')
pause                          


Vector=[];
%Introducing Control Blocks
disp('for i=1:10, ')
disp('    Vector(i)=i; ')
disp('end ')
disp('Vector ')
disp('plot(Vector) ')
for i=1:10,
    Vector(i)=i;
end
Vector
plot(Vector)
shg
pause
disp(' ')
disp(' ')


% NOTE that MATLAB can perfom Matrix manipulation in PARALLEL
% Writing scripts in a "clever" way (in parallel) can speed up
% SIGNIFICANTLY the execution scripts
%

disp('hold on ')
disp('plot((1:10),''r.-'') ')
disp('hold off ')
hold on             % To superimpose different plots, should turn on the "hold" option by typing "hold on"
plot((1:10),'r.-')
hold off
shg
pause
disp(' ')
disp(' ')


disp('Vector=[]; ')
disp('for i=1:10, ')
disp('    if mod(i,2)==0, ') 
disp('        Vector=[Vector, i]; ')
disp('    end ')
disp('end ')
disp('Vector ')
disp('hold on ')
disp('plot((2:2:10),Vector,''ro'') ')
disp('hold off ')
Vector=[];
for i=1:10,
    if mod(i,2)==0, 
        Vector=[Vector, i];
    end
end
Vector
hold on
plot((2:2:10),Vector,'ro')
shg
hold off
pause
disp(' ')
disp(' ')


% Find zero-crossings of a Random Signal with Processing in Parallel
disp('RandSig = randn(1,150); ')
disp('plot(RandSig,''.-'') ')
disp('title(''Random Noise'') ')
RandSig = randn(1,150);
plot(RandSig,'.-')
title('Random Noise')
shg
pause
disp(' ')
disp(' ')

% positive values
disp('hold on ')
disp('Index=find(RandSig(1:end-1)>0); ')
disp('plot(Index,RandSig(Index), ''go'') ')
hold on
Index=find(RandSig(1:end-1)>0);
plot(Index,RandSig(Index), 'go')
pause
disp(' ')
disp(' ')

% zero-crossings
disp('hold on ')
disp('Index=find((RandSig(1:end-1)>0) & (RandSig(2:end)<0)); ')
disp('plot(Index, 0, ''ro'') ')
disp('hold off ')
disp('title(''Random Noise with Zero-crossings'') ')
hold on
Index=find((RandSig(1:end-1)>0) & (RandSig(2:end)<0));
plot(Index, 0, 'ro')
hold off
title('Random Noise with Zero-crossings')
shg
pause
clf
disp(' ')
disp(' ')


%% Introduction To Signal Processing
% 
% Introduction to basic mathematical functions like exp, log and sin/cos
%
% Herein, we introduce basic "Signal and Systems" functions like the auto- (and cross-) correlation tool "xcorr"
% and different windows like hamming and bartlett. Matlab supports a very
% wide variety of such functions. For further information, type "help" 


disp('Dirac=zeros(11,1); ')
disp('Dirac(6)=1; ')
disp('stem((-5:5),Dirac) ')
disp('title(''Discrete-Time Dirac'') ')
disp('xlabel(''Time (in samples)'') ')
Dirac=zeros(11,1);
Dirac(6)=1;
stem((-5:5),Dirac)
title('Discrete-Time Dirac')
xlabel('Time (in samples)')
pause
disp(' ')
disp(' ')


disp('SigLen = 150; ')
disp('t=0:SigLen; ')
SigLen = 150;
t=0:SigLen;
pause
clf
disp(' ')
disp(' ')

% Sinusoid (simple)
% x = sin(2*pi*n/N), Ts=1/N, fo=1Hz 
disp('Signal=sin(2*pi*t/SigLen);')
disp('subplot(2,1,1) ')
disp('plot(t,Signal) ')
disp('subplot(2,1,2) ')
disp('stem(t,Signal,''.'') ')
disp('title(''Discrete-Time Sinusoid'') ')
disp('xlabel(''Time (in samples)'') ')
Signal=sin(2*pi*t/SigLen);
subplot(2,1,1)
plot(t,Signal)
title('Discrete-Time Sinusoid')
subplot(2,1,2)
stem(t,Signal,'.')
xlabel('Time (in samples)')
shg
pause
clf
disp(' ')
disp(' ')


% Parameters of a cosine
clf;
disp('% Parameters of a cosine');
disp('A = 10; f0 = 1000; phi = pi/3;')     % 1. the three parameters, A=amplitude, f0=frequency, and phi=phase angle.
disp('T0 = 1/f0;')                         % 2. the period, T0, is 1 over the frequency.
disp('tt = -2*T0 : T0/40 : 2*T0;')         % 3. tt is the time axis for the plot, start 2 periods before 0 and quit 2 periods after.
disp('Signal = A*cos(2*pi*f0*tt + phi);')  % 4. the values of the cosine are computed.
disp('stem(tt,Signal,''.''); hold on;')      % 5. the plot
disp('plot(tt,Signal,''r''); hold off;')     % 6. the samples
disp('title(''Sinusoid: x(t) = 10 cos(2*pi*1000*t + pi/3)'');') % 7. title
disp('xlabel(''Time (sec)'');')              % 8. label axis 
disp('grid on');       
clf;
A = 10; f0 = 1000; phi = pi/3;     % 1. the three parameters, A=amplitude, f0=frequency, and phi=phase angle.
T0 = 1/f0;                         % 2. the period, T0, is 1 over the frequency.
tt = -2*T0 : T0/40 : 2*T0;         % 3. tt is the time axis for the plot, start 2 periods before 0 and quit 2 periods after.
Signal = A*cos(2*pi*f0*tt + phi);  % 4. the values of the cosine are computed.
stem(tt,Signal,'.'); hold on;      % 5. the plot
plot(tt,Signal,'r'); hold off;     % 6. the samples
title('Sinusoid: x(t) = 10 cos(2*pi*1000*t + pi/3)'); % 7. title
xlabel('Time (sec)');              % 8. label axis 
grid on                            % 9. show the grid
shg
pause
disp(' ')
disp(' ')



% Exponential Signals
disp('Signal1=exp(t/SigLen); ')
disp('Signal2=exp(2*t/SigLen); ')
disp('title(''Plot of Exponential Fuctions'') ')
disp('legend(''exp(t)'',''exp(2t)'') ')
hold on
Signal1=exp(-0.1*t/SigLen);
Signal2=exp(2*t/SigLen);
subplot(2,1,1)
plot(Signal1)
hold on 
plot(Signal2,'r')
hold off
title('Plot of Exponential Fuctions')
legend('exp(t)','exp(2t)')
subplot(2,1,2)
hold on 
stem((1:10:151),Signal1(1:10:end))
stem((1:10:151),Signal2(1:10:end),'r')
shg
pause
clf
disp(' ')
disp(' ')


% Log Signals
disp('Signal1=log(t/SigLen); ')
disp('Signal2=log(2*t/SigLen); ')
disp('legend(''log(t)'',''log(2t)'') ')
disp('hold on ')
disp('subplot(2,1,1) ')
disp('plot(Signal1) ')
disp('plot(Signal2,''r'') ')
disp('hold off ')
disp('subplot(2,1,2) ')
disp('hold on ')
disp('stem((1:10:151),Signal1(1:10:end)) ')
disp('stem((1:10:151),Signal2(1:10:end),''r'') ')
disp('hold off ')
subplot(2,1,1)
Signal1=log(t/SigLen);
plot(Signal1)
hold on
Signal2=log(2*t/SigLen);
plot(Signal2,'r')
hold off
title('Plot of Log-Fuctions')
legend('log(t)','log(2t)')
subplot(2,1,2)
hold on
stem((1:10:151),Signal1(1:10:end))
stem((1:10:151),Signal2(1:10:end),'r')
hold off
shg
pause
disp(' ')
disp(' ')

% Sinc function
disp('t1 = -20:1/2:20;')
disp('Signal1 = sinc(t1);')
disp('Signal2 = sinc(t1/4);')
disp('Signal3 = sinc(t1/8);')
disp('subplot(2,1,1)')
disp('hold on')
disp('plot(t1,Signal3,''g'')')
disp('plot(t1,Signal2,''r'')')
disp('plot(t1,Signal1)')
disp('hold off')
disp('title(''Plot of Sinc-Fuctions'')')
disp('legend(''sinc(t)'',''sinc(t/4)'',''sinc(t/8)'')')
disp('subplot(2,1,2)')
disp('hold on')
disp('stem(t1,Signal3,''g'')')
disp('stem(t1,Signal2,''r'')')
disp('stem(t1,Signal1)')
disp('hold off')
clf
t1 = (-40:1:40)/2;
Signal1 = sinc(t1);
Signal2 = sinc(t1/4);
Signal3 = sinc(t1/8);
subplot(2,1,1)
hold on
plot(t1,Signal3,'g')
plot(t1,Signal2,'r')
plot(t1,Signal1)
hold off
title('Plot of Sinc-Fuctions')
legend('sinc(t)','sinc(t/4)','sinc(t/8)')
subplot(2,1,2)
hold on
stem(t1,Signal3,'g')
stem(t1,Signal2,'r')
stem(t1,Signal1)
hold off
shg
pause
disp(' ')
disp(' ')


% Cosine of Different frequencies
disp('Signal1=cos(2*pi*t/SigLen*2); ')
disp('Signal1=cos(2*pi*t/SigLen*3); ')
disp('Signal1=cos(2*pi*t/SigLen*4); ')
disp('subplot(2,1,1) ')
disp('hold on ')
disp('plot(Signal1,''r'') ')
disp('plot(Signal2) ')
disp('plot(Signal3,''g'') ')
disp('title(''Plot of Cosines with Different Frequencies'') ')
disp('legend(''f=2 Hz'',''f=3 Hz'',''f=4 Hz'') ')
disp('hold off ')
disp('subplot(2,1,2) ')
disp('hold on ')
disp('stem((1:10:151),Signal1(1:10:end)) ')
disp('stem((1:10:151),Signal2(1:10:end),''r'') ')
disp('stem((1:10:151),Signal3(1:10:end),''g'') ')
disp('hold off ')
SigLen = 150;
t=0:SigLen;
Signal1=cos(2*pi*t/SigLen*2);
Signal2=cos(2*pi*t/SigLen*3);
Signal3=cos(2*pi*t/SigLen*4);
clf
subplot(2,1,1)
hold on
plot(Signal1,'r')
plot(Signal2)
plot(Signal3,'g')
title('Plot of Cosines with Different Frequencies')
legend('f=2 Hz','f=3 Hz','f=4 Hz')
hold off
subplot(2,1,2)
hold on
stem((0:5:150),Signal1(1:5:end),'r')
stem((0:5:150),Signal2(1:5:end))
stem((0:5:150),Signal3(1:5:end),'g')
hold off
shg
disp(' ')
%disp(' END HERE FOR TODAY ')
%pause
disp(' ')
%disp(' END HERE FOR TODAY ')



% % Auto-correlation of a Signal 
% clf
% disp('plot(xcorr(Signal1)) ')
% disp('subplot(2,1,1) ')
% disp('plot(Signal1) ')
% disp('title(''Signal1 in Time Domain'') ')
% disp('subplot(2,1,2) ')
% disp('plot(xcorr(Signal1)) ')
% disp('title(''Auto-correlation of the Signal'') ')
% subplot(2,1,1)
% plot(Signal1)
% title('Signal1 in Time Domain')
% subplot(2,1,2)
% plot(xcorr(Signal1))
% title('Auto-correlation of the Signal')
% shg
% pause
% disp(' ')
% disp(' ')
% 
% 
% clf
% disp('plot(xcorr(Signal2)) ')
% disp('subplot(2,1,1) ')
% disp('plot(Signal2) ')
% disp('title(''Signal2 in Time Domain'') ')
% disp('subplot(2,1,2) ')
% disp('plot(xcorr(Signal2),''r'') ')
% disp('title(''Auto-correlation of the Signal'') ')
% subplot(2,1,1)
% plot(Signal2)
% title('Signal2 in Time Domain')
% subplot(2,1,2)
% plot(xcorr(Signal2),'r')
% title('Auto-correlation of the Signal')
% shg
% pause
% disp(' ')
% disp(' ')
% 
% 
% clf
% disp('plot(xcorr(Signal2),''r'') ')
% disp('plot(xcorr(Signal1+Signal2)) ')
% disp('title(''Auto-correlations of Two Signals'') ')
% disp('legend(''Signal2'',''Signal1+Signal2'') ')
% plot(xcorr(Signal2),'r')
% hold on
% plot(xcorr(Signal1+Signal2))
% title('Auto-correlations of Two Signals')
% legend('Signal2','Signal1+Signal2')
% hold off
% shg
% pause
% disp(' ')
% disp(' ')
% 
% % Cross-correlation of a Signal1 and Signal2 
% disp('plot(xcorr(Signal1,Signal2)) ')
% disp('title(''Cross-correlation of Signal1 and Signal2'') ')
% clf
% plot(xcorr(Signal1,Signal2))
% title('Cross-correlation of Signal1 and Signal2')
% shg
pause
disp(' ')
disp(' ')

clf
% A hamming window is chosen
disp('wRect = rectwin(SigLen); ')
disp('plot(wRect)')
disp('title(''Rectangular Window of 150 samples length'') ')
wRect = rectwin(SigLen);
plot(wRect)
title('Rectangular Window of 150 samples length')
shg
pause
disp(' ')
disp(' ')

clf
% A hamming window is chosen
disp('wHamm = hamming(SigLen); ')
disp('plot(wHamm)')
disp('title(''Hamming Window of 150 samples length'') ')
wHamm = hamming(SigLen);
plot(wHamm)
title('Hamming Window of 150 samples length')
shg
pause
disp(' ')
disp(' ')

% A bartlett window is chosen
% similar to a triangular window as returned by the triang function. 
% However, the Bartlett window always has zeros at the first and last samples, 
% while the triangular window is nonzero at those points.
disp('wBart = bartlett(SigLen);')
disp('plot(wBart)')
disp('title(''Bartlett Window of 150 samples length'') ')
wBart = bartlett(SigLen);
plot(wBart)
title('Bartlett Window of 150 samples length')
shg
pause
disp(' ')
disp(' ')


%% SAMPLING
% x(t) = 20cos(2pi(40)t-0.4pi); 40Hz CT sinusoid;

clf
disp('f0 = 40;')
disp('T0 = 1/f0;')
disp('Ts = [5 2.5 0.1]/1000;')
disp('for i=1:3')
disp('    tn = (-T0 : Ts(i) : T0);')  
disp('    xn = 20*cos(2*pi*f0*tn-0.4*pi);')
disp('    subplot(3,1,i);')
disp('    stem(tn,xn,''r.'');hold on;')
disp('    plot(tn,xn);') 
disp('    title([''Samples of Sinusoid: Ts = '' num2str(Ts(i)) '' sec'']);')
disp('end;')
f0 = 40;
T0 = 1/f0;
Ts = [5 2.5 0.1]/1000;
for i=1:3
    tn = (-T0 : Ts(i) : T0);  
    % n = tn/Ts(i);
    xn = 20*cos(2*pi*f0*tn-0.4*pi);
    subplot(3,1,i);
    stem(tn,xn,'r.');hold on;
    plot(tn,xn); 
    title(['Samples of Sinusoid: Ts = ' num2str(Ts(i)) ' sec']);
end;
shg
pause
disp(' ')
disp(' ')


%% FFT BLOCK

close all
disp('playshow fft')
%playshow fftdemo
pause
disp(' ')
disp(' ')


SigLen = 151;
t=1:SigLen;


disp('Signal1=cos(2*pi*t/SigLen*2); ')
disp('Signal3=cos(2*pi*t/SigLen*4); ')
Signal1=cos(2*pi*t/SigLen*2);
Signal3=cos(2*pi*t/SigLen*4);


clf
% fft of a signal: Magnitude and Angle

disp('% DFT of Signal1');
disp('fftSignal=fft(Signal1); ')
disp('subplot(3,1,1) ')
disp('plot(Signal1) ')
disp('title(''Time-Domain Signal'') ')
disp('subplot(3,1,2) ')
disp('plot(20*log(abs(fftSignal))) ')
disp('ylabel(''Magnitude (In db)'') ')
disp('subplot(3,1,3) ')
disp('plot(angle(fftSignal)) ')
disp('ylabel(''Angle'') ')
disp('xlabel(''Frequency'') ')
fftSignal=fft(Signal1);
subplot(3,1,1)
plot(Signal1)
title('Time-Domain Signal')
subplot(3,1,2)
plot(20*log(abs(fftSignal)))
ylabel('Magnitude (In db)')
subplot(3,1,3)
plot(angle(fftSignal))
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause 


clf
% fft of Signal1
disp('% Increase DFT Frequencies (sample points/FFT length)...');
fftSignal=fft(Signal1,512);
subplot(3,1,1)
plot(Signal1)
title('Time-Domain Signal')
subplot(3,1,2)
plot(20*log(abs(fftSignal(2:end))))
ylabel('Magnitude (In db)')
subplot(3,1,3)
plot(angle(fftSignal))
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause 

disp('% Highlight (in red) conjugate symmetry...');
subplot(3,1,2)
hold on
plot(20*log(abs(fftSignal(2:end/2))),'r')
title('DFT of the Signal 1')
ylabel('Magnitude (In db)')
subplot(3,1,3)
hold on
plot(angle(fftSignal(1:end/2)),'r')
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause

clf
disp('% Use fftshift to shift zero-frequency component to center of spectrum ...');
fftSignal=fft(Signal1,512);
fftSignal=fftshift(fftSignal(2:end));

subplot(3,1,1)
plot(Signal1)
title('Time-Domain Signal 1')
subplot(3,1,2)
plot(20*log(abs(fftSignal(2:end))))
ylabel('Magnitude (In db)')
subplot(3,1,3)
plot(angle(fftSignal))
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause 

disp('% Show final Signal Spectrum...');
subplot(3,1,2)
hold on
plot((256:509),20*log(abs(fftSignal(end/2+1:end-1))),'r')
title('DFT of the Signal 1')
subplot(3,1,3)
hold on;
plot((257:510),angle(fftSignal(end/2+1:end-1)),'r')
disp(' ')
disp(' ')
shg
pause


clf
% fft of Signal3
disp('% DFT of Signal3');
fftSignal=fft(Signal3,512);
fftSignal=fftshift(fftSignal(2:end));

subplot(3,1,1)
plot(Signal3)
title('Time-Domain Signal 1')
subplot(3,1,2)
plot(20*log(abs(fftSignal(2:end))))
ylabel('Magnitude (In db)')
subplot(3,1,3)
plot(angle(fftSignal))
ylabel('Angle')
xlabel('Frequency')
subplot(3,1,2)
hold on
plot((256:509),20*log(abs(fftSignal(end/2+1:end-1))),'r')
title('DFT of the Signal 1')
subplot(3,1,3)
hold on;
plot((257:510),angle(fftSignal(end/2+1:end-1)),'r')
disp(' ')
disp(' ')
shg
pause



clf
% fft of Signal1+Signal3
disp('% DFT of Signal1+Signal3');
fftSignal=fft(Signal1+Signal3,512);
subplot(3,1,1)
plot(Signal1+Signal3)
title('Time-Domain Signal 1 + Signal 3')
subplot(3,1,2)
plot(20*log(abs(fftSignal(2:end/2))))
title('DFT of the Signal1+Signal3')
ylabel('Magnitude (In db)')
subplot(3,1,3)
plot(angle(fftSignal(1:end/2)))
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause


clf
% fft of Signal1+Signal3
disp('% DFT of Signal1+Signal3+Random Noise');
fftSignal=fft(Signal1+Signal3+rand(1,SigLen),512);
subplot(3,1,1)
plot(Signal1+Signal3+rand(1,SigLen))
title('Time-Domain Signal Plus Random Noise')
subplot(3,1,2)
plot(20*log(abs(fftSignal(2:end/2))))
ylabel('Magnitude (In db)')
title('DFT of the Signal1+Signal3 Plus Random Noise')
subplot(3,1,3)
plot(angle(fftSignal(1:end/2)))
ylabel('Angle')
xlabel('Frequency')
disp(' ')
disp(' ')
shg
pause


%% Fourier Transform / Gibbs effect
% Changing the number of Discrete Frequencies (Fourier length)
clf
disp('% Changing the number of Discrete Frequencies (Fourier length)');
t1 = (-40:1:40)/2;
SignalS = sinc(t1/2);
N=[128 256 512];
figure
subplot(4,1,1);
stem(t1,SignalS,'r.'); axis tight; hold on;
plot(t1,SignalS);
for i=1:3
    fftSignal=fft(SignalS,N(i));
    fftSignal=fftshift(fftSignal(1:end));
    subplot(4,1,i+1)
    stem(abs(fftSignal),'r.'); axis tight; hold on;
    plot(abs(fftSignal))
    title(['DFT of sinc function: N = ' num2str(N(i)) ' samples']);
end;
disp(' ')
disp(' ')
shg
pause

disp('% Changing the number of Samples (Signal length)');
t1 = (-40:1:40)/2;
N=512;
T = [40 100 400];
figure
for i=1:3
    t1 = (-T(i):1:T(i))/2;
    SignalS = sinc(t1/2);
    fftSignal=fft(SignalS,N);
    fftSignal=fftshift(fftSignal(1:end));
    subplot(3,1,i)
    stem(abs(fftSignal(150:370)),'r.'); axis tight; hold on;
    plot(abs(fftSignal(150:370)))
    title(['DFT of sinc function: T = ' num2str(2*T(i)) ' samples']);
end;



%% Convolution in the Time and Frequency Domains
% 

warning off
SigLen = 151;
t=1:.35:SigLen;


disp(' ')
disp(' ')
disp('Signal1=sin(2*pi*t.^2/SigLen); ')
disp('Signal2=exp(-3.5*(-1:.15:1).^2); ')
disp('ConvSigTime = conv(Signal1,Signal2); ')
Signal1=sin(2*pi*t.^2/SigLen);
Signal2=exp(-3.5*(-1:.15:1).^2);
ConvSigTime = conv(Signal1,Signal2);
disp(' ')
disp(' ')

disp('subplot(3,1,1) ')
disp('plot(Signal1) ')
disp('title(''Signal1'') ')
disp('subplot(3,1,2) ')
disp('plot(Signal2,''r'') ')
disp('title(''Signal2'') ')
disp('subplot(3,1,3) ')
disp('plot(ConvSigTime) ')
disp('title(''Convolution of 2 Signals in the Time Domain'') ')
disp('xlabel(''Time (in Samples)'') ')
disp(' ')
disp(' ')


clf
subplot(3,1,1)
plot(Signal1)
title('Signal1')
shg
pause 
subplot(3,1,2)
plot(Signal2,'r')
title('Signal2')
pause 
subplot(3,1,3)
plot(ConvSigTime)
title('Convolution of 2 Signals in the Time Domain')
xlabel('Time (in Samples)')
shg
pause


% Convolution in the Frequency Domain
% FFT with sufficient number of samples
disp('Convolution in the Frequency Domain')
disp('FFT with Sufficient number of samples')
disp(' ')
disp(' ')
disp('Signal1_FFT=fft(Signal1,512); ')
disp('Signal2_FFT=fft(Signal2,512); ')
disp('ConvSigFreq=Signal1_FFT.*Signal2_FFT; ')
disp('ConvSigTime_Freq=ifft(ConvSigFreq); ')
Signal1_FFT=fft(Signal1,512);
Signal2_FFT=fft(Signal2,512);
ConvSigFreq=Signal1_FFT.*Signal2_FFT;
ConvSigTime_Freq=ifft(ConvSigFreq);
disp(' ')
disp(' ')
pause

disp('subplot(3,1,1) ')
disp('plot(20*log(abs(Signal1_FFT))) ')
disp('ylabel(''Magnitude (in db)'') ')
disp('title(''Magnitude of Signal1''s DFT'') ')
disp('subplot(3,1,2) ')
disp('plot(20*log(abs(Signal2_FFT))) ')
disp('ylabel(''Magnitude (in db)'') ')
disp('title(''Magnitude of Signal2''s DFT'') ')
disp('subplot(3,1,3) ')
disp('plot(20*log(abs(ConvSigFreq))) ')
disp('ylabel(''Magnitude (in db)'') ')
disp('title(''Magnitude of Convolved Signals'' DFT'') ')
disp('subplot(3,1,1) ')
disp('title(''Convolution of 2 Signals in the Frequency Domain'') ')
disp('subplot(3,1,3) ')
disp('xlabel(''Frequency'') ')
disp(' ')
disp(' ')
pause


clf
subplot(3,1,1)
plot(20*log(abs(Signal1_FFT)))
ylabel('Magnitude (in db)')
title('Magnitude of Signal1''s DFT')
pause
subplot(3,1,2)
plot(20*log(abs(Signal2_FFT)))
ylabel('Magnitude (in db)')
title('Magnitude of Signal2''s DFT')
pause
subplot(3,1,3)
hold on
plot(20*log(abs(ConvSigFreq)))
ylabel('Magnitude (in db)')
title('Magnitude of Convolved Signals'' DFT')
pause
subplot(3,1,1)
title('Convolution of 2 Signals in the Frequency Domain')
subplot(3,1,3)
xlabel('Frequency')
pause



% shg
% clf
% plot(ConvSigTime)
% hold on
% plot(ConvSigTime_Freq(1:(length(Signal1)+length(Signal2)-1)),'r:')
% hold off
% title('Comparison of the Convolution in Both Domains')
% legend('Convolution in Time','Convolution in Freq.')
% pause


% 
% % Convolution in the Frequency Domain
% % FFT with insufficient number of samples
% disp('Convolution in the Frequency Domain')
% disp('FFT with insufficient number of samples')
% Signal1_FFT=fft(Signal1,min(length(Signal1),length(Signal2)));
% Signal2_FFT=fft(Signal2,min(length(Signal1),length(Signal2)));
% 
% ConvSigFreq = Signal1_FFT.*Signal2_FFT;
% ConvSigTime_Freq = ifft(ConvSigFreq);
% 
% figure;
% subplot(3,1,1)
% plot(20*log(abs(Signal1_FFT)))
% ylabel('Magnitude (in db)')
% title('Magnitude of Signal1''s DFT')
% subplot(3,1,2)
% plot(20*log(abs(Signal2_FFT)))
% ylabel('Magnitude (in db)')
% title('Magnitude of Signal2''s DFT')
% subplot(3,1,3)
% hold on
% plot(20*log(abs(ConvSigFreq)))
% ylabel('Magnitude (in db)')
% title('Magnitude of Convolved Signals'' DFT')
% subplot(3,1,1)
% title('Convolution of 2 Signals in the Frequency Domain')
% subplot(3,1,3)
% xlabel('Frequency')
% 
% figure;
% plot(ConvSigTime)
% hold on
% plot(ConvSigTime_Freq,'r:')
% hold off
% title('Comparison of the Convolution in Both Domains')
% legend('Convolution in Time','Convolution in Freq.')
% 



% %% Signal Processing Demos of MATLAB
% warning off
% close all
% 
% 
% %specgramdemo
% %pause
% %
% disp('sigdemo1')
% sigdemo1
% disp(' ')
% disp(' ')
% pause
% 
% disp('sigdemo2')
% sigdemo2
% disp(' ')
% disp(' ')
% pause
% 
% % disp('fdatool')
% % fdatool
% % disp(' ')
% % disp(' ')
% % pause
% % 
% % disp('playshow gaussfirdesigndemo')
% % playshow gaussfirdesigndemo
% % disp(' ')
% % disp(' ')
% % pause
% % 
% % disp('psddemo')
% % psddemo
% % disp(' ')
% % disp(' ')
% % pause
% 
% %% SPFirst Demos
% %
% % Download "spfirst" toolbox & demos. Follow installation instractions.  
% % More: http://users.ece.gatech.edu/mcclella/SPFirst/Updates/SPFirstMATLAB.html
% 
% cd('C:\MATLABR2008b\toolbox\spfirst');
% spfirst;
% 
% % Fourier Series
% fseriesdemo;
% 
% % Discrete Convolution
% dconvdemo;
% 
% % Continuous Convolution
% cconvdemo;
% 
% % Sampling-Aliasing
% con2dis;
% 
% % CT filtering
% cltidemo;
% 
% % DT filtering
% dltidemo;