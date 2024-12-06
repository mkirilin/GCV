function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = gcv(U, S, V, x, b, m0, m, allSV)
    % Pre-compute frequently used values
    coeffs_all = U' * b;
    s = coeffs_all.^2;
    x_norm = norm(x);
    b_norm = norm(b)^2;
    
    % Handle singular value cases
    if allSV
        s = s(:);
        cs = cumsum(s, 'reverse');
        k_values = (1:m0);
    else
        s = [0; s];
        cs = cumsum(s);
        k_values = (0:m0);
    end
    
    % Vectorized denominator computation
    denom = (1 - k_values/m).^2;
    
    % Efficient current sums computation
    if allSV
        current_sums = cs(k_values+1)' ./ denom;
    else
        current_sums = (b_norm - cs(k_values+1)') ./ denom;
    end
    
    % Find optimal k values
    [~, idx] = min(current_sums);
    k_gcv = k_values(idx);
    
    % Scale coefficients once
    coeffs_scaled = coeffs_all ./ S;
    
    % Compute solutions efficiently
    if k_gcv == 0
        Xgcv = zeros(size(x));
    else
        Xgcv = V(:,1:k_gcv) * coeffs_scaled(1:k_gcv);
    end
    
    % Compute errors progressively
    errs = zeros(1, m0);
    running_sum = zeros(size(x));
    
    % Use efficient matrix multiplication for progressive errors
    for k = 1:m0
        running_sum = running_sum + V(:,k) * coeffs_scaled(k);
        errs(k) = norm(running_sum - x) / x_norm;
    end
    
    % Find optimal k and compute solution
    [~, k_opt] = min(errs);
    if k_opt == 0
        Xopt = zeros(size(x));
    else
        Xopt = V(:,1:k_opt) * coeffs_scaled(1:k_opt);
    end

    % Compute final errors
    err_gcv = norm(Xgcv - x) / x_norm;
    err_opt = norm(Xopt - x) / x_norm;
    
    fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
    fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end