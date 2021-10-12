function [L, R, LAR] = reduce2(A, numMod)
  % Solve EVP for right- and left eigenvectors
  [V,D]  = eig(A);
  lambda = diag(D);

  % Sort the eigenvalues so that the modes comes with increasing frequency
  [~,indx] = sort(abs(lambda));
  lambda   = lambda(indx);
  V        = V(:,indx);

  T = zeros(size(V));
  k = 1;
  while k <= length(lambda)
    if imag(lambda(k)) ~= 0.0
      T(:,k)   = real(V(:,k));
      T(:,k+1) = imag(V(:,k));
      k = k+2;
    else
      T(:,k) = V(:,k);
      k = k+1;
    end
  end

  nvec = 1:8;
  nsel = sort([nvec*2-1 nvec*2]);

  % Define new state-space
  LAR = T\A*T;

  % Modal reduction
  LAR = LAR(nsel,nsel);

  Tr = T;
  Tl = inv(T);
  R  = Tr(:,nsel);
  L  = Tl(nsel,:);
  L  = L.';
end
