function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = gcv(U, S, V,x,b,m0,m)

  % Precompute squared dot products
  s = (U' * b).^2; % k in [1,...,m0]
  s = [0; s]; % Add 0 to the beginning, size of s is m0+1
  if m0 == m
    s = s(1:end-1); % Remove last element
  end

  % Compute cumulative sums
  cs = cumsum(s);

  % Compute norm of b
  b_norm = norm(b)^2;

  % Define k range
  if m0 == m
    k_values = (0:m0-1);
  else
    k_values = (0:m0);
  end

  % Compute denominators for all k
  denom = (1 - (k_values) / m).^2;

  % Compute current sums for all k
  current_sums = (b_norm - cs(k_values + 1)') ./ denom;

  % Find k that minimizes the current sum
  [~, idx] = min(current_sums);
  k_gcv = k_values(idx);

  % Compute cumulative coefficients
  coeffs_all = (U' * b) ./ diag(S);

  % Compute cumulative solutions X for all k <= m0
  X_cumsum = cumsum(V .* coeffs_all', 2);
  
  % Compute errors for each k
  errs = vecnorm(X_cumsum - x, 2, 1);
  
  % Find k that minimizes the error
  [~, k_opt] = min(errs);
  
  % Compute X using vectorized operations
  if k_gcv == 0
    Xgcv = zeros(size(x));
  else
    Xgcv = X_cumsum(:, k_gcv);
  end

  % Compute optimal Xopt
  if k_opt == 0
    Xopt = zeros(size(x));
  else
    Xopt = X_cumsum(:, k_opt);
  end

  err_gcv = norm(Xgcv - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
  fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end