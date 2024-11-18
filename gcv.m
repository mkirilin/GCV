function [X, Xopt, err_gcv, err_opt] = gcv(U, S, V,x,b,m0,m)

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
  argmin_k = k_values(idx);

  % Compute cumulative coefficients
  coeffs_all = (U' * b) ./ diag(S);

  % Compute cumulative solutions X for all k <= m0
  X_cumsum = cumsum(V .* coeffs_all', 2);
  
  % Compute errors for each k
  errs = vecnorm(X_cumsum - x, 2, 1);
  
  % Find k that minimizes the error
  [~, k_opt] = min(errs);
  
  % Compute X using vectorized operations
  X = X_cumsum(:, argmin_k);

  % Compute optimal Xopt
  Xopt = X_cumsum(:, k_opt);

  err_gcv = norm(X - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, argmin_k);
  fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end