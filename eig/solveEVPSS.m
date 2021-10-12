function eigenSolution = solveEVPSS(clSys, varargin)

  % Solve EvP
  [U, Lambda] = eig(clSys.SS.A);

  % Sort eigen- value and vectors
  lambda   = diag(Lambda);
  %[~, idx] = sort(abs(lambda));
  [~, idx] = sort(abs(imag(lambda)));
  lambda   = lambda(idx);
  U        = U(:, idx);

  for i = 1:size(U, 2)
    U(:,i) = U(:,i)/max(U(:, i));
  end

  eigenSolution.lambdas = lambda;

  w_Hz = abs(imag(lambda))/2/pi;

  % Remove real eigenvalues
  %find(w_Hz > 0.5)
  %w_Hz   = w_Hz(w_Hz > 0.5);
  %fprintf('%10.6f Hz\n', w_Hz(1:16))

  % Collect solution for export
  eigenPairs = cell(1, length(lambda));
  for i = 1:length(lambda)
    eigenPairs{i} = {lambda(i), U(:,i)};
  end

  eigenSolution.eigenPairs = eigenPairs;
  eigenSolution.lambdas    = lambda;

  eigenSolution.gist = sprintf('%14.6f rad/s %14.6f Hz\n', ...
                               [abs(imag(lambda(1:2:80, 1)))'; ...
                                w_Hz(1:2:80, 1)']);

  if nargin > 1
    fprintf('Frequency summary:\n')
    disp(eigenSolution.gist)
  end
end
