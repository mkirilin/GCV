function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = gcv(S,x,b,m0,m,allSV,idx,n)

  %coeffs_all = U' * b;

  b_image = reshape(b, [n, n]);
  coeffs_all_func = dct2(b_image);
  coeffs_all_func = coeffs_all_func(:);
  coeffs_all_func = coeffs_all_func(idx);
  %assert(norm(coeffs_all - coeffs_all_func) < 1e-8,...
  %      'Error: coeffs_all and coeffs_all_func are not equal');
  
  % Precompute squared dot products
  s = coeffs_all_func.^2; % k in [1,...,m0]
  if allSV
    s = s(:);
  else
    s = [0; s]; % Add 0 to the beginning, size of s is m0+1
  end
  %if m0 == m
  %  s = s(1:end-1); % Remove last element
  %end

  % Compute cumulative sums
  if allSV
    cs = cumsum(s, 'reverse');
  else
    cs = cumsum(s);
  end

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
  if allSV
    current_sums = cs(k_values + 1)' ./ denom;
  else
    current_sums = (b_norm - cs(k_values + 1)') ./ denom;
  end

  % Find k that minimizes the current sum
  [~, idx_min] = min(current_sums);
  k_gcv = k_values(idx_min);

  % Compute cumulative coefficients
  coeffs_all_func = coeffs_all_func ./ diag(S);

  % Compute cumulative solutions X for all k <= m0
  %X_cumsum = cumsum(V .* coeffs_all', 2);

  % Initialize result
  X_cumsum_func = zeros(n*n, n*n);

  % Get DCT matrix once
  D = dctmtx(n)';
  
  % For each coefficient in sorted order
  for k = 1:n*n
      % Get pre-sorting position
      [i, j] = ind2sub([n,n], idx(k));
      
      % Get basis vectors for original position
      ui = D(:,j);
      uj = D(:,i);
      
      % Create 2D basis vector via kronecker product
      v = kron(ui, uj);
      
      % Multiply by sorted coefficient and store
      X_cumsum_func(:,k) = v * coeffs_all_func(k);
  end
  
  % Compute cumulative sum
  X_cumsum_func = cumsum(X_cumsum_func, 2);

  %assert(norm(X_cumsum(1:1000) - X_cumsum_func(1:1000)) < 1e-8,...
  %      'Error: X_cumsum and X_cumsum_func are not equal');
  
  % Compute errors for each k
  errs = vecnorm(X_cumsum_func - x, 2, 1);
  
  % Find k that minimizes the error
  [~, k_opt] = min(errs);
  
  % Compute X using vectorized operations
  if k_gcv == 0
    Xgcv = zeros(size(x));
  else
    Xgcv = X_cumsum_func(:, k_gcv);
  end

  % Compute optimal Xopt
  if k_opt == 0
    Xopt = zeros(size(x));
  else
    Xopt = X_cumsum_func(:, k_opt);
  end

  err_gcv = norm(Xgcv - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
  fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end