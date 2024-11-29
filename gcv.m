function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] =...
  gcv(S, x, b, m0, m, allSV, idx, n)

  %coeffs_all = U' * b;

  b_image = reshape(b, [n, n]);
  coeffs = dct2(b_image);
  coeffs = coeffs(:);
  coeffs = coeffs(idx);
  %assert(norm(coeffs_all - coeffs_all_func) < 1e-8,...
  %      'Error: coeffs_all and coeffs_all_func are not equal');
  
  % Precompute squared dot products
  s = coeffs.^2; % k in [1,...,m0]
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
  coeffs = coeffs ./ S;

  % Compute cumulative solutions X for all k <= m0
  %X_cumsum = cumsum(V .* coeffs_all', 2);

  % Get DCT matrix once
  D = dctmtx(n)';

  % Pre-compute all basis vectors efficiently
  [I, J] = ind2sub([n,n], idx);
  U = D(:,J);  % n x n^2
  V = D(:,I);  % n x n^2

  % Compute outer products using matrix operations
  X_all = zeros(n^2, n^2);
  for k = 1:n^2
      % Get kth columns
      uk = U(:,k);  % n x 1
      vk = V(:,k);  % n x 1
      
      % Compute outer product
      X_all(:,k) = kron(uk, vk);
  end

  % Instead of computing full cumsum, calculate errors progressively
  errs = compute_chunked_errors(X_all, coeffs, x, n);
  
  % Find optimal k
  [~, k_opt] = min(errs);
  
  % Compute only needed solutions
  if k_gcv == 0
      Xgcv = zeros(size(x));
  else
      Xgcv = X_all(:,1:k_gcv) * coeffs(1:k_gcv);
  end
  
  if k_opt == 0
      Xopt = zeros(size(x));
  else
      Xopt = X_all(:,1:k_opt) * coeffs(1:k_opt);
  end
  
  % Compute errors
  err_gcv = norm(Xgcv - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);

  % Print results
  fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
  fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end

function errs = compute_chunked_errors(V, coeffs, x, n)
  % Initialize with optimal chunk size
  chunk_size = 10;
  errs = zeros(1, n^2);
  running_sum = sparse(n^2, 1);
  x_norm = norm(x);
  
  % Pre-compute chunk boundaries
  chunk_starts = 1:chunk_size:n^2;
  chunk_ends = min(chunk_starts + chunk_size - 1, n^2);
  
  % Process each chunk more efficiently
  for i = 1:length(chunk_starts)
    idx_range = chunk_starts(i):chunk_ends(i);
    
    % Pre-compute V*coeffs for chunk
    chunk_V = V(:,idx_range);
    chunk_coeffs = coeffs(idx_range);
    
    % Process coefficients progressively but more efficiently
    for j = 1:length(idx_range)
      k = idx_range(j);
      running_sum = running_sum + chunk_V(:,j) * chunk_coeffs(j);
      
      % More efficient norm computation using dot product
      diff = running_sum - x;
      errs(k) = sqrt(dot(diff,diff)) / x_norm;
    end
  end
end