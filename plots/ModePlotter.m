classdef ModePlotter

properties
  numDof;    % number of degress of freedom
  numNo;     % number of nodes
  Omega;     % angular velocity [rad/s]
  eigenData; % cell array containing eigen- values and vectors
end


methods
  function obj = ModePlotter(eigenSolution, numDof)
    % Constructor.
    obj.Omega     = eigenSolution.Omega;
    obj.eigenData = eigenSolution.eigenPairs;
    obj.numDof    = numDof;
    obj.numNo     = obj.numDof/4;
  end

  function plotModeShape(obj, modeNumber)
    % Plot a mode from the eigenSolution object.

    % Actual mode (forward and absolute)
    actualMode = 2*modeNumber;

    % Natural frequency and mode-shape of the chosen mode
    wn = imag(obj.eigenData{actualMode}{1});
    U  = obj.eigenData{actualMode}{2};

    % Total simulation time
    t = linspace(0, 8/abs(wn), 100);

    % Select v and w DOFs from the mode shape
    v = U(1:4:obj.numDof);
    w = U(2:4:obj.numDof);

    % Calculate the modal displacements v and w for each node (u)
    u = ndgrid(1:obj.numNo, 1:length(t));
    v = real(v)*cos(wn*t) + imag(v)*sin(wn*t);
    w = real(w)*cos(wn*t) + imag(w)*sin(wn*t);

    figure(); set(gcf,'color','w'); hold on;
    for node = 1:obj.numNo
      plot3(u(node,:), v(node,:), w(node,:), 'k', 'LineWidth', 2.5);
    end

    view(-25, 20); grid on;
    xticks(1:1:obj.numNo)
    xlabel('Node Number'); ylabel('v'); zlabel('w');
    set(findall(gcf, '-property', 'FontSize'), ...
        'FontSize', 15, 'FontName', 'Times')
    title(sprintf('\\Omega: %.2f rpm    Mode: %d    Nat. Freq.: %.2f Hz', ...
          [obj.Omega*60/2/pi, modeNumber, abs(wn/2/pi)]), 'FontSize', 20)
  end
end
end
