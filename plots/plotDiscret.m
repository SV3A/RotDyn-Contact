function plotDiscret(mesh, varargin)
  % Plots the discretization of the shaft.

  showNodes = false;

  if nargin > 1
    if strcmp(varargin{1}, 'show-nodes')
      showNodes = true;
    end
  end

  l  = 1e3 * mesh(1, :);
  rx = 1e3 * mesh(2, :);

  e_start = 0;
  e_end   = 0;

  figure('Name','Discretization');

  for i = 1:length(l)
    e_end = e_end + l(i);

    plot([e_start e_end e_end e_start e_start], ...
         [-rx(i) -rx(i) rx(i) rx(i) -rx(i)], '-k', 'LineWidth', 1.5); hold on
    e_start = e_end;
  end

  grid on
  box off

  xlim([0 sum(l)])
  ylim([-100 100])

  xticks(50:50:sum(l) )
  yticks(-60:15:60 )

  % Plot node numbers
  if showNodes
    for i = 1:length(l)+1
      if mod(i,2) ~= 0 || i == length(l)+1
        text(sum(l(1:i-1)), 0, num2str(i), 'FontSize', 16)
      end
    end
  end

  % Plot elements mounted on the shaft
  %if length(elems) >= 1
  %  for j = 1:length(elems)
  %    scatter(sum(l(1:elems(j)-1)), -rx(elems(j)-1), 260,'or','filled')

  %    text(sum(l(1:elems(j)-1)) - 15, -(rx(elems(j)-1) + 15), el_names{j},...
  %      'FontSize', 26, 'FontWeight', 'bold')
  %  end
  %end

  set(gcf,'color','w');
  set(gca, 'FontUnits', 'normalized', 'FontSize', 0.07, 'FontName', 'Times')
  xlabel('Shaft length [mm]', 'interpreter', 'latex', ...
         'FontUnits', 'normalized', 'FontSize', 0.1)
  ylabel('Radius [mm]', 'interpreter', 'latex', 'FontUnits', 'normalized', ...
         'FontSize', 0.1)
end
