function eigenSolution = solveEVP(rs, Omega, EVPType, varargin)
  % Handle for solving different types of EVPs.
  %
  % Input:
  %   rs:        rotor-system struct exported from the 'RotorFEModel' class.
  %   Omega:     angular velocity [rad/s]
  %   EVPType:   which eigenproblem formulation to use (see switch below for
  %              options).
  %   show-gist: optionally one can supply the trailing string argument:
  %              'show-gist', this shows a small summary of found frequencies.
  %
  % Output:
  %   eigenSolution: cell array containing the solution to the EVP, each cell
  %                  contains a eigenvalue and the associated eigenvector.

  % Call the appropriate analysis
  switch EVPType
    case 'std'
      eigenSolution = stdFormulation(rs, Omega);
    case 'general'
      eigenSolution = genFormulation1(rs, Omega);
    case 'general2'
      eigenSolution = genFormulation2(rs, Omega);
    case 'general3'
      eigenSolution = genFormulation3(rs, Omega);
    case 'stdgen'
      eigenSolution = stdGenFormulation(rs, Omega);
    case 'norot'
      if Omega > 0
        warning('No-rot requested with Omega > 0, ignoring this value.')
      end
      eigenSolution = norotFormulation(rs);
    otherwise
      warning('Unknown EVP-type, analysis failed.')
  end

  eigenSolution.EVPType = EVPType;

  % Show summary of frequencies
  if nargin > 3
    if strcmp(varargin{1}, 'show-gist')
      fprintf('Frequency summary (%s formulation):\n', eigenSolution.EVPType)
      disp(eigenSolution.gist)
    end
  end
end


function eigenSolution = stdFormulation(rs, Omega)
  % Solve the standard EVP: A*phi_i = lambda_i*phi_i

  % Build states matrix A as:
  %     |    0       I    |
  % A = |-M^-1*K  -M^-1*C |
  if isfield(rs, 'D')
    A = [ zeros(rs.numDof)        eye(rs.numDof)
            -rs.M\rs.K     -rs.M\(-Omega*rs.G+rs.D) ];
  else
    A = [ zeros(rs.numDof)        eye(rs.numDof)
            -rs.M\rs.K        -rs.M\(-Omega*rs.G)   ];
  end

  eigenSolution = solveStateSpaceEVP(A, Omega, rs.numDof*2);
end


function eigenSolution = stdGenFormulation(rs, Omega)
  % Solve the standard EVP: A*phi_i = lambda_i*I*phi_i

  % Build states matrix A as:
  %     |    0       I    |
  % A = |-M^-1*K  -M^-1*C | B = I
  if isfield(rs, 'D')
    A = [ zeros(rs.numDof)        eye(rs.numDof)
            -rs.M\rs.K     -rs.M\(-Omega*rs.G+rs.D) ];
  else
    A = [ zeros(rs.numDof)        eye(rs.numDof)
            -rs.M\rs.K        -rs.M\(-Omega*rs.G)   ];
  end

  eigenSolution = solveStateSpaceEVP(A, eye(2*rs.numDof), Omega, rs.numDof*2);
end


function eigenSolution = genFormulation1(rs, Omega)
  % Solves the generalized EVP: -A*phi_i = lambda_i*B*phi_i,
  %
  % with the state-space defined as:
  %   | M  0 ||q_dd|   | C  K ||q_d|   |0|
  %   | 0  K || q_d| + |-K  0 || q | = |F|,  with  [C] = -Omega*[G] + [D]
  % or:
  %      [B] * {z_d} +    [A] * {z}  = {f}
  %
  % inserting the homgeneous solution {z} = {phi} e^(lambda*t) we get:
  %     lambda*[B]{phi} + [A]{phi} = {0}
  %     ==>
  %     -[A]*{phi} = lambda*[B]*{phi}

  zeroM = zeros(rs.numDof);

  % Build state matrices A and B as:
  %      |-C -K |      | M  0 |
  % -A = | K  0 |, B = | 0  K |
  B    = sparse([  rs.M       zeroM
                  zeroM        rs.K ]);

  % Add damping matrix if present in the 'rs' object
  if isfield(rs, 'D')
    negA = sparse([ Omega*rs.G-rs.D -rs.K
                          rs.K      zeroM ]);
  else
    negA = sparse([ Omega*rs.G      -rs.K
                       rs.K         zeroM ]);
  end

  eigenSolution = solveStateSpaceEVP(negA, B, Omega, rs.numDof*2);
end


function eigenSolution = genFormulation2(rs, Omega)
  % Solves the generalized EVP: -A*phi_i = lambda_i*B*phi_i,
  %
  % with the state-space defined as:
  %   | 0  M ||q_dd|   |-M  0 ||q_d|   |0|
  %   | M  C || q_d| + | 0  K || q | = |F|,  with  [C] = -Omega*[G] + [D]
  % or:
  %      [B] * {z_d} +    [A] * {z}  = {f}

  zeroM = zeros(rs.numDof);

  % Build state matrices A and B as:
  %      | M  0 |      | 0  M |
  % -A = | 0 -K |, B = | M  C |
  negA = sparse([  rs.M  zeroM
                  zeroM  -rs.K ]);

  if isfield(rs, 'D')
    B = sparse([ zeroM         rs.M
                  rs.M -Omega*rs.G+rs.D ]);
  else
    B = sparse([ zeroM         rs.M
                  rs.M     -Omega*rs.G    ]);
  end

  eigenSolution = solveStateSpaceEVP(negA, B, Omega, rs.numDof*2);
end


function eigenSolution = genFormulation3(rs, Omega)
  % Solves the generalized EVP: -A*phi_i = lambda_i*B*phi_i,
  %
  % with the state-space defined as:
  %   | C  M || q_d|   | K  0 || q |   |0|
  %   | M  0 ||q_dd| + | 0 -M ||q_d| = |F|,  with  [C] = -Omega*[G] + [D]
  % or:
  %      [B] * {z_d} +    [A] * {z}  = {f}

  zeroM = zeros(rs.numDof);

  % Build state matrices A and B as:
  %      |-K  0 |      | C  M |
  % -A = | 0  M |, B = | M  0 |
  negA = sparse([ -rs.K  zeroM
                  zeroM   rs.M ]);

  if isfield(rs, 'D')
    B = sparse([ -Omega*rs.G+rs.D rs.M
                       rs.M       zeroM ]);
  else
    B = sparse([ -Omega*rs.G      rs.M
                     rs.M         zeroM ]);
  end

  eigenSolution = solveStateSpaceEVP(negA, B, Omega, rs.numDof*2);
end


function eigenSolution = norotFormulation(rs)
  % Solve the non-rotating EVP: [K]{U} = w^2*[M]{U}

  % DEBUG: Print the conditioning number of the state matrices
  %fprintf('Condition numbers:\n  cond M: %14.6f cond K: %14.6f\n\n', ...
  %        [cond(rs.M) cond(rs.K)])

  [U, Lambda] = eig(rs.K, rs.M);

  [lam, sortIdx] = sort(diag(Lambda));
  U = U(:, sortIdx);

  for i = 1:size(U, 2)
    U(:,i) = U(:,i) ./ max(U(:,i));
  end

  w    = sqrt(lam);
  w_Hz = w/2/pi;

  % Collect solution for export
  for i = 1:length(lam)
    eigenPairs{i} = {lam(i), U(:,i)};
  end
  eigenSolution.Omega      = 0;
  eigenSolution.eigenPairs = eigenPairs;

  eigenSolution.gist = sprintf('%14.6f rad/s %14.6f Hz\n', ...
                               [abs(w(1:10, 1))'; w_Hz(1:10, 1)']);
end


function eigenSolution = solveStateSpaceEVP(A, varargin)
  % Common function between the state-space formulations for solving the EVP
  % given the state matrices A or A and B.

  % Determine whether standard or generalized EVP is to be solved
  if nargin < 4
    Omega = varargin{1};
    %nEigs = varargin{2};

    % 'eigs' returns the eigenvalues 'Lambda' and the eigenvectors 'U':
    %   'Lambda' is a 2*numDof x 2*numDof matrix, with the eigenvalues stored
    %   in its diagonal.
    %   'U' is a 2*numDof x 2*numDof matrix
    %[U, Lambda] = eigs(A, nEigs, 'sm');
    [U, Lambda] = eig(A);
  else
    B     = varargin{1};
    Omega = varargin{2};
    nEigs = varargin{3};

    [U, Lambda] = eigs(A, B, nEigs, 'sm');
  end

  % DEBUG: Print the conditioning number of the state matrices
  %if exist('B')
  %  fprintf('Conditions numbers:\n  cond -A: %14.6f cond B: %14.6f\n\n', ...
  %          [condest(A) condest(B)])
  %else
  %  fprintf('Conditions number:\n  cond -A: %14.6f\n\n', condest(A))
  %end

  % Sort eigen- value and vectors
  lambda   = diag(Lambda);
  %[~, idx] = sort(abs(lambda));
  [~, idx] = sort(abs(imag(lambda)));
  lambda   = lambda(idx);
  U        = U(:, idx);

  % Normalize the eigenvectors
  %for i = 1:size(U, 2)
  %  U(:,i) = U(:,i) ./ max(U(:,i));
  %end

  % Natural frequency
  w_Hz = abs(imag(lambda))/2/pi;

  % Collect solution for export
  eigenPairs = cell(1, length(lambda));
  for i = 1:length(lambda)
    eigenPairs{i} = {lambda(i), U(:,i)};
  end
  eigenSolution.Omega      = Omega;
  eigenSolution.eigenPairs = eigenPairs;
  eigenSolution.lambdas    = lambda;

  %lambda = lambda(w_Hz > 0.5);
  %w_Hz   = w_Hz(w_Hz > 0.5);

  eigenSolution.gist = sprintf('%14.6f rad/s %14.6f Hz\n', ...
                               [abs(imag(lambda(1:2:20, 1)))'; ...
                                w_Hz(1:2:20, 1)']);
end
