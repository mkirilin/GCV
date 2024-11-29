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

  % Multiply by coefficients and compute cumulative sum
  X_cumsum = cumsum(X_all .* coeffs', 2);

  %assert(norm(X_cumsum(1:1000) - X_cumsum_func(1:1000)) < 1e-8,...
  %      'Error: X_cumsum and X_cumsum_func are not equal');
  
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

function [V] = compute_V(n, idx)
  % Compute V for computing results by the formula:
  % X = cumsum(V .* coeffs', 2);

  % Get DCT matrix once
  D = dctmtx(n)';

  % Pre-compute all basis vectors efficiently
  [I, J] = ind2sub([n,n], idx);
  U = D(:,J);  % n x n^2
  V = D(:,I);  % n x n^2

  % Pre-allocate arrays for sparse matrix construction
  nz_per_col = n;  % Each column will have n non-zero elements
  total_nz = n^2 * nz_per_col;
  
  % Create index arrays for sparse matrix
  rows = zeros(total_nz, 1);
  cols = zeros(total_nz, 1);
  vals = zeros(total_nz, 1);
  
  % Fill index arrays efficiently
  idx = 1;
  for k = 1:n^2
    uk = U(:,k);
    vk = V(:,k);
    
    % Compute outer product indices
    [r, c] = ndgrid(1:n, k);
    curr_vals = uk .* vk;
    
    % Store in arrays
    num_elements = numel(r);
    rows(idx:idx+num_elements-1) = r(:);
    cols(idx:idx+num_elements-1) = c(:);
    vals(idx:idx+num_elements-1) = curr_vals(:);
    
    idx = idx + num_elements;
  end
  
  % Construct sparse matrix directly
  V = sparse(rows, cols, vals, n^2, n^2);
end

function [k_opt, err_gcv, err_opt, Xgcv, Xopt] =...
  compute_errors(V, coeffs, x, k_gcv, n)

  % Compute cumulative solutions X for all k <= m0 by the formula:
  % X = cumsum(V .* coeffs', 2);

  % Instead of computing full cumsum, calculate errors progressively
  errs = zeros(1, n^2);
  running_sum = sparse(n^2, 1);

  % Progressive reconstruction and error calculation
  for k = 1:n^2
    running_sum = running_sum + V(:,k) * coeffs(k);
    errs(k) = norm(running_sum - x);
  end

  % Find optimal k
  [~, k_opt] = min(errs);

  % Compute only needed solutions
  if k_gcv == 0
    Xgcv = zeros(size(x));
  else
    Xgcv = V(:,1:k_gcv) * coeffs(1:k_gcv);
  end

  if k_opt == 0
    Xopt = zeros(size(x));
  else
    Xopt = V(:,1:k_opt) * coeffs(1:k_opt);
  end

  % Compute errors
  err_gcv = norm(Xgcv - x) / norm(x);
  err_opt = norm(Xopt - x) / norm(x);
end

function errs = compute_chunked_errors(V, coeffs, x, n)
  % Initialize
  chunk_size = min(1000, n^2);  % Adjust chunk size based on available memory
  num_chunks = ceil(n^2/chunk_size);
  errs = zeros(1, n^2);
  running_sum = sparse(n^2, 1);
  
  % Process chunks
  for i = 1:num_chunks
      % Calculate chunk indices
      start_idx = (i-1)*chunk_size + 1;
      end_idx = min(i*chunk_size, n^2);
      chunk_range = start_idx:end_idx;
      
      % Compute chunk sum efficiently
      chunk_sum = V(:,chunk_range) * coeffs(chunk_range);
      running_sum = running_sum + chunk_sum;
      
      % Compute errors for chunk
      errs(chunk_range) = sqrt(sum((running_sum - x).^2));
  end
end