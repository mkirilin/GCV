function [Xhdp, err_hdp, k_hdp] =...
  hdpBlurGauss(S, x, bn, m0, m, allSV, idx, n)

  bn_image = reshape(bn, [n, n]);
  coeffs = dct2(bn_image);
  coeffs = coeffs(:);
  coeffs = coeffs(idx);

  % Precompute squared dot products
  s = coeffs.^2; % k in [1,...,m0]
  if allSV
    s = s(:);
  else
    s = [0; s]; % Add 0 to the beginning, size of s is m0+1
  end

  % Compute cumulative sums
  if allSV
    cs = cumsum(s, 'reverse');
  else
    cs = cumsum(s);
  end

  % Define k range
  if m0 == m
    k_values = (0:m0-1);
  else
    k_values = (0:m0);
  end

  % Compute denominators for all k
  denom = (S(1:m0).^2)';

  % Compute current sums for all k
  if allSV
    current_sums = cs(k_values + 1)' ./ denom;
  else
    current_sums = (norm(bn)^2 - cs(k_values + 1)') ./ denom;
  end

  % Find k that minimizes the current sum
  [~, idx_min] = min(current_sums);
  k_hdp = k_values(idx_min);

  % Compute cumulative coefficients
  coeffs = coeffs ./ S;

  % Compute cumulative solutions X for all k <= m0
  %X_cumsum = cumsum(V .* coeffs_all', 2);

  % Get DCT matrix once
  D = dctmtx(n)';

  % Pre-compute all basis vectors efficiently
  [I, J] = ind2sub([n,n], idx);
  U = D(:,I);  % n x n^2
  V = D(:,J);  % n x n^2

  % Set chunk size based on available memory
  chunk_size = 128;  % Adjust as needed

  % Compute solutions for k_gcv
  Xhdp = compute_solution_vectorized(U, V, coeffs, k_hdp, n, chunk_size);

  % Compute errors
  err_hdp = norm(Xhdp - x) / norm(x);

  % Deprecated. Uncomment if needed.
  % Print results
  % fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
end

function X_rec = compute_solution_vectorized(U, V, coeffs, k_max, n, chunk_size)
  if k_max == 0
    X_rec = zeros(n^2, 1);
    return;
  end

  total_k = k_max;
  num_chunks = ceil(total_k / chunk_size);
  X_rec = zeros(n^2, 1);

  for chunk = 1:num_chunks
    start_k = (chunk - 1) * chunk_size + 1;
    end_k = min(chunk * chunk_size, total_k);
    indices = start_k:end_k;
    chunk_len = length(indices);

    % Extract chunks
    U_chunk = U(:, indices);           % n x chunk_len
    V_chunk = V(:, indices);           % n x chunk_len
    coeffs_chunk = coeffs(indices);    % chunk_len x 1

    % Compute outer products vectorized
    U_reshaped = reshape(U_chunk, [n, 1, chunk_len]);   % n x 1 x chunk_len
    V_reshaped = reshape(V_chunk, [1, n, chunk_len]);   % 1 x n x chunk_len

    % Compute outer products
    outer_products = bsxfun(@times, U_reshaped, V_reshaped);  % n x n x chunk_len
    outer_products = reshape(outer_products, [n^2, chunk_len]);

    % Multiply by coefficients and sum
    X_rec = X_rec + outer_products * coeffs_chunk;
  end
end