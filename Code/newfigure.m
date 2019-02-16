function fig = newfigure(figname,xinch,yinch)
fig = figure;
set(gcf,'Name',figname);
set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf, 'PaperSize', [xinch yinch]);
set(gcf,'Position',[3,4,xinch,yinch]);
set(gcf,'PaperUnits','Inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperPosition',get(gcf,'Position'));