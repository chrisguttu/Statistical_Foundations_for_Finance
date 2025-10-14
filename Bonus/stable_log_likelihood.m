%% creating the log-likelihood
function logL = stable_log_likelihood(a, b, x)
    % Ensure parameters are within valid range
    if a <= 0 || a >= 2 || abs(b) > 1
        error('Invalid parameters for stable distribution');
    end

    % Calculate the log-PDF for each data point in x using asymstabpdf
    pdf_values = asymstabpdf(x, a, b);
    log_pdf_values = log(pdf_values);

    % Compute the log-likelihood as the sum of log-PDF values
    logL = sum(log_pdf_values);
end
