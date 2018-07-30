%Experiment code to analyse data from the images. Pre-processed using gaia.
%Adam Dad


clc;
clear all;
close all;

M = csvread('allstars1.csv'); %reads from the file allstars.csv

MJD = M(1,:); %creates vectors of all the numbers from M
MJD = MJD-MJD(1); % subtract fist element from array from the rest gives days

for j = 2:2631 %iteration of code for all stars
data = M(j,:); 

DateBins = linspace(MJD(1),MJD(98),51); % groups the data in even x-axis groups

%Takes the mean value in each group of data produces the lines fitted to scatter:
MovingMean = zeros(length(DateBins)-1); 
  for i = 1:length(DateBins)-1
      low = DateBins(i);
      high = DateBins(i+1);
      mask = (MJD>low) & (MJD<high);
      InBin = data(mask);
      MeanBin = mean(InBin);
      MovingMean(i) = MeanBin;
  end

centers = MJD(1)+DateBins(2)-DateBins(1):DateBins(2)-DateBins(1):MJD(98); % center of each bin
%{
figure; %commented out to prevent graphs of all data being produced
scatter(MJD,data,'.')
hold on;
plot(centers, MovingMean)
hold off; 
%}
Fs = 1/(DateBins(2)-DateBins(1));

Y = fft(data); %fast furiour transform of the data
P2 = abs(Y/(length(DateBins)-1)); 
P1 = P2(1:length(centers)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(length(centers)/2))/length(centers);% Above lines required to produce the fourier transform.
period = 1./f;
%figure;
%plot(period,P1) 
mask = period < 15; 
pks = findpeaks(P1(mask));
index = P1 == max(pks);
P(j,:) = period(index);
Av(j,:) = mean(data);
stdev(j,:) = std(data) ;
end

% Okay now for the variables and tranients:

figure;
scatter(Av, stdev, '.') %produces a scatter of mean against standard deviation, allowing transients to be seen.
stdevT = find(stdev > 0.7);

%following code produces graphs for a specifec set:
for k = 2:length(stdevT) %iteration of code for data set
data = M(k,:);

DateBins = linspace(MJD(1),MJD(98),51);

MovingMean = zeros(length(DateBins)-1);
  for i = 1:length(DateBins)-1
      low = DateBins(i);
      high = DateBins(i+1);
      mask = (MJD>low) & (MJD<high);
      InBin = data(mask);
      MeanBin = mean(InBin);
      MovingMean(i) = MeanBin;
  end

centers = MJD(1)+DateBins(2)-DateBins(1):DateBins(2)-DateBins(1):MJD(98);

figure;
scatter(MJD,data,'.')
hold on;
plot(centers, MovingMean)
hold off;
Fs = 1/(DateBins(2)-DateBins(1));

Y = fft(data);
P2 = abs(Y/(length(DateBins)-1));
P1 = P2(1:length(centers)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(length(centers)/2))/length(centers);
period = 1./f;
%following plot produces fourier transforms for the periods
%figure;
%plot(period,P1) 
mask = period < 15;
pks = findpeaks(P1(mask));
index = P1 == max(pks);
P(k,:) = period(index);
%Av = mean(data) 
%stdev= std(data)
end

