pathSetup

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


rotMod.addNodeComponent(6, disc)

rotMod.addNodeComponent(9, pmbMass)
rotMod.addNodeComponent(9, pmbStiff, 'internal')
rotMod.addNodeComponent(9, pmbDamper)

rotMod.addNodeComponent(12, magTarget)

rotMod.addNodeComponent(15, sphBearingStiff, 'internal')
rotMod.addNodeComponent(15, sphBearingDamp)

rotMod.addNodeComponent(17, coupling)


rotMod.printInfo()

