function [T_B, T_D] = mkTransMatrices(numDof, discs, bearings)

  numBearDof = 2*length(bearings);

  % Transformation matrix from bearing DOFs to global dofs
  T_B = zeros(numDof, numBearDof);
  for i = 1:length(bearings)
    es = (bearings{i}.nodalPosition - 1)*4;

    offset = (i-1)*2; % equals: 0, 2, 4, 6 ...
    T_B(es+1, offset+1) = 1;
    T_B(es+2, offset+2) = 1;
  end

  % Transformation matrix selecting the discs lateral DOFs
  T_D = zeros(numDof, 2);
  es  = (discs{1}.nodalPosition - 1)*4;
  T_D(es+1, 1) = 1;
  T_D(es+2, 2) = 1;
end
