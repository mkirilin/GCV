function [X] = IRgcv(A,b)
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
  k_values = (0:floor(m/2))';

  % Compute denominators for all k
  denom = (1 - k_values / m).^2;

  % Compute current sums for all k
  current_sums = cs_rev(k_values + 1) ./ denom;

  % Find k that minimizes the current sum
  [~, idx] = min(current_sums);
  argmin_k = k_values(idx);

  disp(argmin_k);

  % Compute coefficients
  coeffs = (U(:, 1:argmin_k)' * b) ./ diag(S(1:argmin_k, 1:argmin_k));

  % Compute X using vectorized operations
  X = V(:, 1:argmin_k) * coeffs;
end