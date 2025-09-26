function [Xhdp, err_hdp, k_hdp] = ...
  hdp(S, V, x, bn, m0, m, allSV, coeffs_all)

  s = coeffs_all.^2;
  coeffs_scaled = coeffs_all ./ S(1:m0+1);

  x_norm = norm(x);
  k_values = (0:m0);

  denom = (S(1:m0+1).^2)';
  fprintf('Size of denom: %d x %d\n', size(denom, 1), size(denom, 2));

  % Compute current sums
  if (allSV)
    current_sums = compute_current_sums_allSV(s, k_values, denom);
  else
    current_sums = compute_current_sums(s, bn, k_values, denom);
  end

  % Compute ALL ERRORS progressively
  errs = zeros(1, m0);
  running_sum = zeros(size(x));
  for k = 1:m0
    running_sum = running_sum + V(:,k) * coeffs_scaled(k);
    errs(k) = norm(running_sum - x) / x_norm;
  end

  % Find gcv k value and compute gcv solution
  [~, idx] = min(current_sums);
  k_hdp = k_values(idx);
  Xhdp = V(:,1:k_hdp) * coeffs_scaled(1:k_hdp);
  if k_hdp == 0 err_hdp = 1; else err_hdp = errs(k_hdp); end

  % Deprecated. Uncomment if needed.
  %fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
  %fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
end

function current_sums = compute_current_sums_allSV(s, k_values, denom)
  s = s(:);
  cs = cumsum(s, 'reverse');
  cs = [cs; 0];

  current_sums = cs(k_values+1)' ./ denom;
end

function current_sums = compute_current_sums(s, bn, k_values, denom)
  s = s(:);
  cs = cumsum(s);
  cs = [0; cs];

  current_sums = (norm(bn)^2 - cs(k_values+1)') ./ denom;
end