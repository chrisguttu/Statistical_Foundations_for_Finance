%Modified Hartley
function [solvec, crit, iter] = hartley_mod(y, init, tol, maxit)
    k = size(init, 1) / 3; % Assuming init has 3 rows per component: mean, std, and weight
    old = init;
    new = zeros(size(init));
    iter = 0;
    crit = 0;
    
    while 1
        iter = iter + 1;
        
        % Extracting parameters for each component
        mu = old(1:k);
        sigma = old(k+1:2*k);
        lambda = old(2*k+1:end);
        
        % E-step: Calculating responsibilities
        gamma = zeros(length(y), k);
        for i = 1:k
            gamma(:, i) = lambda(i) * normpdf(y, mu(i), sigma(i));
        end
        gamma = gamma ./ sum(gamma, 2);
        
        % M-step: Updating parameters
        Nk = sum(gamma, 1);
        new(1:k) = sum(y .* gamma) ./ Nk; % Updating means
        for i = 1:k
            new(k+i) = sqrt(sum(gamma(:, i) .* (y - new(i)).^2) / Nk(i)); % Updating std devs
        end
        new(2*k+1:end) = Nk / length(y); % Updating weights
        
        % Convergence check
        crit = max(abs(old - new));
        solvec = new;
        if any(isnan(solvec)), break, end
        if (crit < tol) || (iter >= maxit), break, end
        
        old = new;
    end
end
