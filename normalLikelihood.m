function [ logLikelihood ] = normalLikelihood( params,data )
%normalLikelihood - returns the (inverted) likelihood of observing 'data' if drawing
%from a normal distribution with parameters mu and sigma.

mu = params(1);
sigma = params(2);

logLikelihood = -1.*sum(log(normpdf(data,mu,sigma)));% Inverted so that minimization algorithms can function properly
end
