%This code is provided to implement the Michaud(1998) statistical
%resampling technique for optimal portfolio construction.  

%The data has been provided, which includes all stocks which were in the S&P 500 from
%April 2001 (when stocks were first quoted to the penny) to March 2007 (the
%latest data available on the CRSP website as of June 2009); however, if you would like
%to update the data, this code will handle that as well (data can be
%obtained by students at wrds.wharton.upenn.edu).  It is important,
%however, that the data take the form [permno,date,return,return+1,log(return)] 
%in order to work properly.  

clear all;
clc;
close all;

RawData=xlsread('MichaudDataLogReturns');

%Discover how many months data we are dealing with, where Nn will be the first
%occurance of a new permno:
FirstPermno=RawData(1,1);
Nn=1;
while RawData(Nn,1)==FirstPermno;
    Nn=Nn+1;
end

%n is the number of months
n=Nn-1;

%How many stocks do we have?
M=length(RawData(:,1))/n;

%Now, since the princomp command needed to perform Principal Component
%Analysis (PCA) requires us to have each column to represent a variable (stock) and
%each row to represent an observaion (month), we must set up our data
%matrix accordingly.  

TimeSeries=zeros(n,M);
for i=1:M;
    for k=1:n;
        TimeSeries(k,i)=RawData((n*(i-1))+k,3);
    end
end

%Annualize the returns:
TimeSeries=TimeSeries.*12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now that we have our data set up, we perform a quick PCA to draw down our
%data to as few components as we wish; we will try 10 components.  Note
%that this exactly mimicks the code written by Jason Hsu for these
%purposes.

NumComp=10;
[eigenvect,prin_score,eigenval]=princomp(TimeSeries);
PCAData=zeros(n,NumComp);
for i=1:NumComp;
    [loadings, positions] = sort(eigenvect(:,i),'descend');
    long_assets = sort(positions(1:round(M/10)));
    short_assets = sort(positions(round(M*(9/10)):end));
    long_returns = mean(TimeSeries(:,long_assets),2);
    short_returns = mean(TimeSeries(:,short_assets),2);
    PCAData(:,i) = long_returns - short_returns;
end

%From now on, we will be working with this data instead of our original
%time series.  PCAData is composed of return sequences for the 10 factor
%mimicking portfolios.

SampleMean=mean(PCAData);
SampleSigma=cov(PCAData);

%We may find it useful to run our initial MV optimization once here using
%the sampled moments.  Also, note that NumMVF is the number of points we
%will calculate along each MVF we construct.

NumMVF=50;
[OriginalWeight,OriginalMean,OriginalSD]=...
    mean_var_frontier(SampleMean',SampleSigma,NumMVF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We are finally ready to begin generating our resampled return sequences.
%We will set NumResampled to be the number of times we would like to resample and
%recalculate the MVF.

NumResampled=20;
ResampledData=zeros(n,NumComp,NumResampled);
ResampledMean=zeros(1,NumComp,NumResampled);
ResampledSigma=zeros(NumComp,NumComp,NumResampled);
ResampledMVFWeight=zeros(NumMVF,NumComp,NumResampled);
ResampledMVFMean=zeros(1,NumMVF,NumResampled);
ResampledMVFSD=zeros(1,NumMVF,NumResampled);

for i=1:NumResampled;
    ResampledData(:,:,i)=mvnrnd(SampleMean,SampleSigma,n);
    ResampledMean(:,:,i)=mean(ResampledData(:,:,i));
    ResampledSigma(:,:,i)=cov(ResampledData(:,:,i));
    [ResampledMVFWeight(:,:,i),ThrowAway1,ThrowAway2]=...
        mean_var_frontier(ResampledMean(:,:,i)',ResampledSigma(:,:,i),NumMVF);
    ResampledRealizedData=PCAData*ResampledMVFWeight(:,:,i)';
    ResampledMVFMean(:,:,i)=mean(ResampledRealizedData);
    ResampledMVFSD(:,:,i)=sqrt(var(ResampledRealizedData));
end

%We will now find our resampled optimal portfolio weights, which are the
%resampled weights averaged across each resample draw.  That is, if we
%resample and construct the MVF 3 times, then our optimal resampled
%portfolio is the one given by averaging the weights on the first
%component, second component, etc. until we have weights for all of the
%components.  These are to be thought of as an average of NumResampled
%statistically equivalent portfolios.

OptimalResampledWeight=mean(ResampledMVFWeight,3);
OptimalResampledRealizedData=PCAData*OptimalResampledWeight';
OptimalResampledMean=mean(OptimalResampledRealizedData);
OptimalResampledSD=sqrt(var(OptimalResampledRealizedData));

%Plot the Original MVF in solid red, the optimal resampled MVF in blue, 
%and the resampled MVFs in dashed yellow.

figure(1)
plot(OriginalSD,OriginalMean,'color','r')
axis([OriginalSD(1,1),OriginalSD(NumMVF,1),OriginalMean(1,1),OriginalMean(NumMVF,1)])
xlabel('Standard Deviation')
ylabel('Mean')
title('Markowitz MVF')
legend('MVF','location','Best')

figure(2)
plot(OriginalSD,OriginalMean,'color','r')
line(OptimalResampledSD,OptimalResampledMean)
for i=1:NumResampled;
    line(ResampledMVFSD(:,:,i),ResampledMVFMean(:,:,i),'LineStyle','--','color','y')
end
axis([OriginalSD(1,1),OriginalSD(NumMVF,1),OriginalMean(1,1),OriginalMean(NumMVF,1)])
xlabel('Standard Deviation')
ylabel('Mean')
title('MVF and Statistically Equivalent Resampled MVFs')
legend('Original MVF','MRE MVF','Statistically Equivalent MVFs','location','Best')





