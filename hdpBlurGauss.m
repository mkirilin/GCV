function [k_hdp] =...
  hdpBlurGauss(S, bn, m0, idx, n)

  bn_image = reshape(bn, [n, n]);
  coeffs = dct2(bn_image);
  coeffs = coeffs(:);
  coeffs = coeffs(idx);

  % Take only first m0 singular values
  coeffs = coeffs(1:m0);

  % Precompute squared dot products
  s = coeffs.^2; % k in [1,...,m0]
  s = s(:);

  % Compute cumulative sums
  cs = cumsum(s, 'reverse');
  cs = [cs; 0];

  % Define k range
  k_values = (0:m0-1);

  % Compute denominators for all k
  denom = (S(1:m0).^2)';

  % Compute current sums for all k
  current_sums = cs(k_values + 1)' ./ denom;

  % Find k that minimizes the current sum
  [~, idx_min] = min(current_sums);
  k_hdp = k_values(idx_min);

end
