function y = mixnormsim(mu, sig, lam, n, seed)
    if nargin < 5, seed = 0; end
    k = length(mu);
    if k <= 1, error('At least 2 components'), end
    if abs(sum(lam) - 1) > 1e-10, error('Lambda does not sum to one'), end
    if any(lam <= 0), error('Lambda out of range'), end
    if any(lam >= 1), error('Lambda out of range'), end
    if any(sig <= 0), error('Sigma out of range'), end

    pick = zeros(n, k);
    if seed ~= 0, rng(seed, 'twister'); end % Updated for modern MATLAB versions
    for i = 1:n
        mult = randmultinomial(lam);
        pick(i, mult) = 1;
    end

    X = [];
    for j = 1:k
        if seed ~= 0
            rng(j * seed + 1, 'twister');
            Z = norminv(rand(n, 1));
        else
            Z = randn(n, 1);
        end
        X = [X, mu(j) + sig(j) * Z];
    end

    X = pick .* X;
    y = sum(X, 2);
end

function mult = randmultinomial(lam)
    % Assumes lam sums to 1 and lam is a row vector
    cumlam = [0, cumsum(lam)];
    r = rand();
    mult = find(r > cumlam, 1, 'last');
end
