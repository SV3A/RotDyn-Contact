function figHandle = plotDisc(t, y, outputConfig)
  % Plot simulations of wheel simulations in the global system.
  %
  % Input:
  %   t:             Time vector.
  %   y:             Simulations output vector (length(t) X numOutputs).
  %   outputConfig:  String cell-array from SS.OutputName.
  %
  % Ouput:
  %   figHand:       Figure handle.


  % Labels denoting the displacement of the wheel - to be matched with the
  % output configuration of the state-space model.
  labels = {'z_Dx', 'z_Dy'};

  % Match and pair labels with the output configuration, i.e. signals will be
  % plotted according to the labels list.
  [~, ia, ib] = intersect(outputConfig, labels);
  outIdxs = ia(ib);

  % Convert meter to micro meter
  y(:, outIdxs) = y(:, outIdxs)*1e6;

  % Determine y limits
  %maxDisp =  ceil(max(y(:, outIdxs), [], 'all'));
  %minDisp = floor(min(y(:, outIdxs), [], 'all'));

  figHandle = figure('name', 'DW');
  for i = 1:2
    j = i*2-1; % Odd idx: 1, 3, 5 ...

    subplot(1, 4, j:j+1) % Subplot is 1x4 instead of 1x2 to improve margins
    plot(t, y(:, outIdxs(i)), '-k', 'LineWidth', 1.5); grid on

    %if maxDisp > minDisp;  ylim([minDisp maxDisp]);  end

    title(['Wheel displacement in ' labels{i}(4) '-direction'])

    if i == 1
      ylabel('Displacemnt [$\mu$m]', 'interpreter', 'latex')
    end

    xlabel('Time [s]')
  end

  % Set background color, font, and fontsize
  defPltSettings(gcf);
end
