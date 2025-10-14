%%% TESTING THE CODE

% Generate a bivariate non-central Student's t dataset
rng(42);  % Set seed for reproducibility

% Parameters
nobs = 1000;
df = 5;  % degrees of freedom
mu = [1.5, -0.5];  % mean vector
Sigma = [2, 0.8; 0.8, 1];  % covariance matrix
ncp = [0.5, 1];  % non-centrality parameters

% Generate data
data = mvtrnd(Sigma, df, nobs) + mu + ncp;

% Test the function
[param, stderr, iters, loglik, Varcov] = MVTestimation(data);
disp('Estimated Parameters:');
disp(param);
disp('Standard Errors:');
disp(stderr);
disp(['Log-Likelihood: ', num2str(loglik)]);
disp(['Number of Iterations: ', num2str(iters)]);
