function [L, R, LAR] = reduce(A, numMod)
  % Solve EVP for right- and left eigenvectors
  [R, Lambda, L] = eig(A);

  L = conj(L);

  % Sort eigen- value and both vectors
  lambda   = diag(Lambda);
  [~, idx] = sort(abs(lambda));
  lambda   = lambda(idx);
  R        = R(:, idx);
  L        = L(:, idx);

  % Truncate modes
  R = R(:, 1:numMod);
  L = L(:, 1:numMod);

  % Normalize the eigenvectors
  for i = 1:size(R, 2)
    if L(:,i).' * R(:,i) == 0
      error('Invalid eigenvectors')
    elseif L(:,i).' * R(:,i) < 0
      L(:,i) = -L(:,i);
    end
    scale = sqrt(L(:,i).' * R(:,i));
    L(:,i) = L(:,i) ./ scale;
    R(:,i) = R(:,i) ./ scale;
  end

  %for i = 1: size(R, 2)
  %  L(:,i).'*R(:,i)
  %end

  % Truncate to nm x nm
  LAR = L.'*A*R;
  %LAR = eye(length(LAR)) .* diag(LAR);
end
