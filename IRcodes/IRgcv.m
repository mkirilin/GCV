function [X, Xopt] = IRgcv(A,x,b,m0)
  % Calculate execution time of this function
  tic;
  % Find singular vectors of sparce matrix A
  m = size(A,1);
  assert(m0 <= m, 'm0 must be less than or equal to the number of rows of A');
  fprintf('m = %d, m0 = %d\n', m, m0);
  [U,S,V] = svds(A, m0);
  fprintf('SVDS Execution time: %f\n', toc);

  % Precompute squared dot products
  s = (U' * b).^2;

  % Compute cumulative sums in reverse order
  cs_rev = cumsum(s);

  % Compute norm of b
  b_norm = norm(b)^2;

  % Define k range
  k_values = (0:m0-1);

  % Compute denominators for all k
  denom = (1 - k_values / m).^2;

  % Compute current sums for all k
  current_sums = (b_norm - cs_rev(k_values + 1)') ./ denom;

  % Find k that minimizes the current sum
  [~, idx] = min(current_sums);
  argmin_k = k_values(idx);

  % Compute cumulative coefficients
  coeffs_all = (U' * b) ./ diag(S);

  % Compute cumulative solutions X for all k
  X_cumsum = cumsum(V .* coeffs_all', 2);
  
  % Compute errors for each k
  errs = vecnorm(X_cumsum - x, 2, 1);
  
  % Find k that minimizes the error
  [err_opt, k_opt] = min(errs);
  
  % Compute X using vectorized operations
  X = X_cumsum(:, argmin_k);

  % Compute optimal Xopt
  Xopt = X_cumsum(:, k_opt);

  err_gcv = norm(X - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, argmin_k);
  fprintf('Optimal error: %f on k = %d\n', err_opt, k_opt);
  fprintf('Execution time: %f\n', toc);
end