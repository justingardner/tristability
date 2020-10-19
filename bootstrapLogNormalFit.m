function [mu,sigma,logLikelihood] = bootstrapLogNormalFit(data,titleStr,nBootstraps,plotThings)
% This function takes a [1,n] vector of data to be fit by a log normal
% distribution. It performs nBootstraps fits to bootstrapped samples
% (n randomly selected samples from the data with replacement) and returns
% the parameters found for each fit. A figure shows the mean parameters on
% top of the data, and a second figure shows the mu and sigma parameters found
% for each bootstrap.

data = data(data~=0); %Remove 0's which are considered invalid
data = log(data); %Does log transform on the data
fits = zeros(3,nBootstraps); %Prepare to store fitted parameters mu, sigma and likelihoods
bootstrapIndices = randi(size(data,2),nBootstraps,size(data,2)); %Generate indexes to randomly sample from data
bootstrapSamples = data(bootstrapIndices(:,:)); %Retrieve random samples from data

%For each bootstrapped sample data set, find the best-fitting a and b parameters based on the
%resulting log likelihood of the set.
%options = optimset('Display', 'off');
for i=1:nBootstraps
    [fits(1:2,i), fits(3,i)] = fminsearch(@normalLikelihood,[2 2],[],bootstrapSamples(i,:));
end

%Extract the means and standard deviations of the fitted parameters
mu = mean(fits(1,:));
stdMu = std(fits(1,:));
sigma = mean(fits(2,:));
stdSigma = std(fits(2,:));
logLikelihood = -1.*mean(fits(3,:));

%Plot the data (with normalized counts so that the fitted distribution will
%match)
[counts, binValues] = histcounts(data,ceil(4*sqrt(size(data,2))));
normalizedCounts = counts / (sum(counts)*(binValues(2)-binValues(1)));
if(plotThings == true)
    figure();
    bar(binValues(2:end), normalizedCounts, 'barwidth', 1);
end
X=0:.1:60;%max(data);

%Plot the fitted distribution over the data
pdf = normpdf(X,mu,sigma);    
if(plotThings == true)
    hold on;
    plot(X,pdf,'r','lineWidth',2);
    title(['Best normal fit for ',titleStr,' mu = ',num2str(mu),' stdMu = ',num2str(stdMu),' sigma = ',num2str(sigma),' stdSigma = ',num2str(stdSigma)],'FontSize',16)
    set(gca,'FontSize',16)
    xlim([0 60]);
    xlabel('Duration (seconds)');
    ylabel('Proportion of total events');
end

%Plot all of the a and b parameters found.
if(plotThings == true)
    figure();
    plot(fits(1,:),'r.');
    hold on;
    plot(fits(2,:),'b.');
    legend('Mu params','Sigma params');
end
end
