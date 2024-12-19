function [Xgcv, Xopt, err_gcv, err_opt, k_gcv, k_opt] = gcv(U, S, V, x, b, m0, m, allSV)
  % Pre-compute frequently used values
  coeffs_all = U' * b;
  s = coeffs_all.^2;
  x_norm = norm(x);
  b_norm_2 = norm(b)^2;

  s = s(:);
  cs = cumsum(s);
  cs = [0; cs];

  fprintf('cs(m0+1) - ||b||^2 = %f\n', abs(cs(m0+1) - b_norm_2));

  k_values = (0:m0);

  denom = (1 - k_values/m).^2;

  current_sums = (b_norm_2 - cs(k_values+1)') ./ denom;
  
  % Find optimal k values
  [~, idx] = min(current_sums);
  k_gcv = k_values(idx);

  % Scale coefficients once
  coeffs_scaled = coeffs_all ./ S;

  % Compute ALL ERRORS progressively
  errs = zeros(1, m0);
  running_sum = zeros(size(x));

  for k = 1:m0
    running_sum = running_sum + V(:,k) * coeffs_scaled(k);
    errs(k) = norm(running_sum - x) / x_norm;
  end
  
  % Find optimal k and compute solution
  [~, k_opt] = min(errs);
  if (k_opt == 0)
    Xopt = zeros(size(x));
  else
    Xopt = V(:,1:k_opt) * coeffs_scaled(1:k_opt);
  end

  % Compute GCV solution
  if (k_gcv == 0)
    Xgcv = zeros(size(x));
  else
    Xgcv = V(:,1:k_gcv) * coeffs_scaled(1:k_gcv);
  end

  % Compute final errors
  err_gcv = norm(Xgcv - x) / x_norm;
  err_opt = norm(Xopt - x) / x_norm;

  fprintf('GCV error: %f on k = %d\n', err_gcv, k_gcv);
  fprintf('Optimal error: %f on k = %d\n\n', err_opt, k_opt);
  if (allSV)
    compute_errors_allSV(s, V, k_values, denom, x, x_norm, coeffs_scaled);
  end
end

function compute_errors_allSV(s, V, k_values, denom, x, x_norm, coeffs_scaled)
  s_all = s(:);
  cs_all = cumsum(s_all, 'reverse');
  cs_all = [cs_all; 0];

  current_sums_all = cs_all(k_values+1)' ./ denom;

  [~, idx_all] = min(current_sums_all);
  k_gcv_all = k_values(idx_all);

  if (k_gcv_all == 0)
    Xgcv_all = zeros(size(x));
  else
    Xgcv_all = V(:,1:k_gcv_all) * coeffs_scaled(1:k_gcv_all);
  end

  err_gcv_all = norm(Xgcv_all - x) / x_norm;
  fprintf('GCV error (all): %f on k = %d\n', err_gcv_all, k_gcv_all);
end