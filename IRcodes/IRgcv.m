function [X] = IRgcv(A,x,b)
  % Find singular vectors of sparce matrix A
  m = size(A,1);
  disp(m);
  [U,S,V] = svds(A, m, "largest");

  % TODO: take leading singular values k <= m0 = #singular

  % Precompute squared dot products
  s = (U' * b).^2;

  % Compute cumulative sums in reverse order
  cs_rev = cumsum(s, 'reverse');

  % Define k range
  k_values = (0:floor(m/2));

  % Compute denominators for all k
  denom = (1 - k_values / m).^2;

  % Compute current sums for all k
  current_sums = cs_rev(k_values + 1)' ./ denom;

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

  err_gcv = norm(X - x);
  err_opt = norm(Xopt - x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, argmin_k);
  fprintf('Optimal error: %f on k = %d\n', err_opt, k_opt);
end