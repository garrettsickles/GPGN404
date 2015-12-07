clear all;
close all;

%% Test of synthetic (for making sure your code is doing what is expected) 
nt = 200;
t = sort(0+0.5*rand(nt,1));
xt = 2*cos(10*pi*t+pi/6) + cos(20*pi*t+pi/3); %Generate some smooth data
xt = (0.8+0.4*rand(nt,1)).*xt; % Add a lot of noise to it

nn = 51;
Ts = 0.5/(nn-1);
tn = 0.0:Ts:0.5;

xn = xnfromxt(xt,t,tn,Ts,1e-1);
figure; plot(t,xt,'r.',tn,xn,'b')


%% Temperature data %%
data = importdata('Lab6_t_T.csv');
data = data.data;
dates = data(:, 1);
TmaxOrig = data(:, 2);
TminOrig = data(:, 3);
nOrig = size(TmaxOrig);

figure; plot(1:nOrig,TmaxOrig,'r',1:nOrig,TminOrig,'b')

% Tmin
inData = [transpose(0:nOrig-1) TminOrig];

% Do interpolation for some chunk of data
i1=1;
i2=nOrig;
inData1 = inData(i1:i2,:);
inData1(inData1(:,2)==-9999.0,:) = [];

Ts = 4;

xn = xnfromxt(inData1(:,2),inData1(:,1),i1-1:Ts:i2-1,Ts,1e-1);

figure; plot(i1-1:i2-1,inData(i1:i2,2),'r.',i1-1:Ts:i2-1,xn,'b')

