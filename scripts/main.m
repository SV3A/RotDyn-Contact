clear all
close all

%/ Setup /
pathSetup

% Handle Octave specific setup
global isOctave = exist ('OCTAVE_VERSION', 'builtin');

if isOctave
  pkg load control
end


%/ Analyses flags /
fSimulate   = 1;
fReduce     = 1;
fEigVals    = 1;
fLambdas    = 1;
fModeShapes = 1;


%/ Script Setup /
Omega = 17*2*pi/60; % [rad/s]
unbalance = 0.0767*0.001; % [kg/m]

% Define shaft mesh and material

% Shaft discretization:
%
%     Disc
%      __
%     |  |     __   Mag. target
%  ___|  |____|  |______---_______-
% |___|  |____|  |______   _______ |
%     |  |    |__|      ---       -
%     |__|                     Coupling
%           PMB rotor
%

% Length [mm]
% Outer radius [mm]
% Inner radius [mm]
% Partition num


shaftDim = [50.0 119.0 88.0 99.2 72.5 63.3
            12.5  12.5 12.5 12.5 12.5 12.5
             0.0   0.0  0.0  0.0  0.0  0.0
               2     3    3    3    3    2];


msh = Mesh(shaftDim);

% Visualize discretization
%plotDiscret(msh.elements, 'show-nodes');

msh.setDensity(2600);
msh.setEmod(71e9);

% Initiate FE model
rotMod = RotorFEModel(msh.elements);
rotMod.addRayDamping(0, 2.4795e-6);


% Define machine elements
disc = Disc(0.250, 219584.55e-9, 426774.25e-9, unbalance);

pmbMass  = Disc(0.560, 287431.88e-9, 521378.74e-9, 0);

pmbStiff = Bearing([3.09e4    0
                        0   3.09e4]);

if Omega > 0
  pmbDamper = Damper(8.48);

else
  pmbDamper = Damper(40.9);

end

magTarget  = Disc(0.2566, 79854.29e-9, 63014.72e-9, 0);

sphBearingStiff = Bearing([1e9   0
                            0   1e9]);
sphBearingDamp = Damper(100);

coupling = Disc(0.429, 172380e-9, 240578.44e-9, 0);

% FREE-FREE overwrite
%pmbStiff.localK = [1e3   0
%                    0   1e3];
%sphBearingStiff.localK = [1e3   0
%                           0   1e3];


rotMod.addNodeComponent(6, disc)

rotMod.addNodeComponent(9, pmbMass)
rotMod.addNodeComponent(9, pmbStiff, 'internal')
rotMod.addNodeComponent(9, pmbDamper)

rotMod.addNodeComponent(12, magTarget)

rotMod.addNodeComponent(15, sphBearingStiff, 'internal')
rotMod.addNodeComponent(15, sphBearingDamp)

rotMod.addNodeComponent(17, coupling)


rotMod.printInfo()

% Export rotor and clean up
rotSys = rotMod.export();
delete(rotMod);


%/ Build state-space model /

% Construct transformation matrices
[~, T_D] = mkTransMatrices(rotSys.numDof, rotSys.discs, rotSys.bearings);

if isfield(rotSys, 'D');  C = -Omega*rotSys.G+rotSys.D;
else;                     C = -Omega*rotSys.G     ; end

% System matrix
A = [ zeros(rotSys.numDof) eye(rotSys.numDof)
       -rotSys.M\rotSys.K    -rotSys.M\C      ];

% Input: unbalance force
B = [ zeros(rotSys.numDof, 2)
          rotSys.M\T_D        ];

% Output: Displacement and velocity at bearings, and displacement at disc
C = [ T_D.' zeros(2, rotSys.numDof) ];

% Feedforward
D = zeros(size(C, 1), size(B, 2));


if fReduce == 1
  [L, R, A] = reduce2(A, 16);
  rotSS = ss(A, L.'*B, C*R, D);

  for i = 1:size(rotSS.StateName, 1)
    rotSS.StateName{i} = ['eta_' num2str(i)];
  end

  redSys.L = L;
  redSys.R = R;

else
  rotSS = ss(A, B, C, D);


  for i = 1:rotSys.numDof
    rotSS.StateName{i} = ['q_' num2str(i)];
    rotSS.StateName{i+rotSys.numDof} = ['qd_' num2str(i)];
  end
end


rotSS.InputName  = {'funb_x' 'funb_y'};
rotSS.OutputName = {'z_Dx' 'z_Dy'};


if fSimulate == 1
  % Simulation parameters
  dt    = 1e-6;
  tspan = [0 1];

  redSys.SS = rotSS;

  % Initial conditions
  % Displacemnt/velocity and current
  q0 = zeros(2*rotSys.numDof, 1);

  % State vector
  x0 = [redSys.L.' * q0];

  t = (tspan(1):dt:tspan(2))';

  % Unbalance input
  unbx = disc.u*Omega^2*cos(Omega*t);
  unby = disc.u*Omega^2*sin(Omega*t);
  u = [unbx unby];

  [y, ~, x] = lsim(redSys.SS, u, t, x0);

  %[P1, fftFrq] = mkFFT(1/dt, size(y, 1), y(:,1)); plotFFT(P1*1e6, fftFrq);

  pltWDis = plotDisc(t, y, redSys.SS.OutputName);
end


% Perform modal analysis

if fEigVals || fLambdas || fModeShapes
  es1 = solveEVP(rotSys, Omega, 'general', 'show-gist');
  es2 = solveEVP(rotSys, Omega, 'std',     'show-gist');

  if fLambdas == 1
    es1.lambdas
    es2.lambdas
  end

  if fModeShapes == 1
    modPlt = ModePlotter(es1, rotSys.numDof);

    modPlt.plotModeShape(2);
    modPlt.plotModeShape(3);
    modPlt.plotModeShape(6);
    modPlt.plotModeShape(8);
  end
end

