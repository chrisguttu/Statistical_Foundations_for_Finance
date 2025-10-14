% EM for k-mixture of normals. 
function [parameters, logLikelihood, iterations] = mixnormEM(y, k, tol, maxit)
    if nargin < 4, maxit = 5e4; end
    if nargin < 3, tol = 1e-6; end

    % Initialization
    n = length(y);
    mu = linspace(min(y), max(y), k);
    sigma = repmat(std(y) / sqrt(k), 1, k);
    lambda = ones(1, k) / k;
    parameters = [mu; sigma; lambda];
    
    logLikelihood = -inf;
    iterations = 0;
    
    % EM algorithm
    while iterations < maxit
        iterations = iterations + 1;
        
        % E-step
        gamma = zeros(n, k);
        for i = 1:k
            gamma(:, i) = lambda(i) * normpdf(y, mu(i), sigma(i));
        end
        gamma = gamma ./ sum(gamma, 2);
        
        % M-step
        Nk = sum(gamma, 1);
        mu = sum(y .* gamma) ./ Nk;
        sigma = sqrt(sum(gamma .* (y - mu).^2) ./ Nk);
        lambda = Nk / n;
        
        % Log likelihood
        newLogLikelihood = sum(log(sum(gamma, 2)));
        
        % Check for convergence
        if abs(newLogLikelihood - logLikelihood) < tol
            break
        end
        logLikelihood = newLogLikelihood;
        
        % Update parameters
        parameters = [mu; sigma; lambda];
    end
end
