function defPltSettings(figHandle, fontSize)
  % Utility that sets background color, font, and font-size for plots.
  % The font-size can be controller globally by setting the integer variable
  % 'GLOBAL_PLT_FONTSIZE' in the base workspace.
  %
  % The function either takes an optional figure handle or defaults to the
  % current active figure via the 'gcf' object.

  % Set bg color to white
  set(gcf, 'color', 'w');

  % Check if fontSize supplied, or global variable is set, or default to 20,
  % then apply font-size to all text elements in the supplied figure handle or
  % the current figure.
  if ~exist('fontSize', 'var')
    try
      fontSize = evalin('base', 'GLOBAL_PLT_FONTSIZE');
    catch
      fontSize = 20;
    end
  end

  if ~exist('figHandle', 'var')
    figHandle = gcf;
  end

  set(findall(figHandle, '-property', 'FontSize'), 'FontSize', fontSize, ...
      'FontName', 'Times')
end
