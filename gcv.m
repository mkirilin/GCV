function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = ...
  gcv(S, V, x, bn, m0, m, allSV, coeffs_all)

  s = coeffs_all.^2;
  coeffs_scaled = coeffs_all ./ S;

  x_norm = norm(x);
  k_values = (0:m0);
  denom = (1 - k_values/m).^2;

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
  k_gcv = k_values(idx);
  Xgcv = V(:,1:k_gcv) * coeffs_scaled(1:k_gcv);
  if k_gcv == 0 err_gcv = 1; else err_gcv = errs(k_gcv); end

  % Find optimal k and compute optimal solution
  [~, k_opt] = min(errs);
  Xopt = V(:,1:k_opt) * coeffs_scaled(1:k_opt);
  if k_opt == 0 err_opt = 1; else err_opt = errs(k_opt); end

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