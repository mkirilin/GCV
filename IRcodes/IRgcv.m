function [X] = IRgcv(A,b)
% Cast sparse matrix to dense matrix
  Afull = full(A);

% Find singular vectors of sparce matrix A
  [U,S,V] = svd(Afull);

  m = size(Afull,1);
  disp(m);
  min_sum = Inf;
  argmin_k = 0;

  for k = 0:m/2
    dot_products = 0;
    for j = k+1:m
      dot_products = dot_products + (U(:,j)'*b)^2;
    end
    current_sum = dot_products / (1 - k/m)^2;
    if current_sum < min_sum
      min_sum = current_sum;
      argmin_k = k;
    end
  end

  disp(argmin_k);

  X = 0;
  for j = 1:argmin_k
    X = X + (U(:,j)'*b)/(S(j,j))*V(:,j);
  end
end