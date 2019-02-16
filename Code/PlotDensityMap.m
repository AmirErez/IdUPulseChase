function densityout = PlotDensityMap(expdata,tp,ytick,xtick,toPlotPCA)

densityout = struct;
densityout.angl = NaN;
densityout.mean = NaN;
densityout.median = NaN;
densityout.density = expdata.CoarseGrid.ScatterSliceCellCounts{tp};

imagesc(expdata.CoarseGrid.x, expdata.CoarseGrid.y,...
    log10(densityout.density));

hold on
set(gca,'YDir','Normal');
%             colormap(bone(256));
colormap(jet(256));

set(gca,'FontSize',7);
%             t=title(char(expdata.Type(tp)));
%             t.Interpreter = 'None';
%             set(gcf,'Name',['ScatterSlice Density ' char(expdata.Type(tp))]);

%             ytick = expdata.CoarseGrid.y(1:2:end);

set(gca,'YTick',ytick);
yticklab = cell(length(ytick),1);
for yt=1:length(yticklab)
    yticklab{yt} = ['10^{' num2str(ytick(yt)) '}'];
end
set(gca,'YTickLabel',yticklab);

set(gca,'XTick',xtick);
xticktab = cell(length(xtick),1);
for xt=1:length(xticktab)
    xticktab{xt} = ['10^{' num2str(xtick(xt)) '}'];
end
set(gca,'XTickLabel',xticktab);
xlim([min(xtick) max(xtick)]);
%             set(gca,'XTick',expdata.CoarseGrid.x(1:4:end));
%                     set(gca,'XTickLabel',cellstr(['10^{' num2str(expdata.CoarseGrid.x(1:4:end)) '}']));
%
ylim([min(ytick) max(ytick)]);
% plot(expdata.CoarseGrid.x( expdata.CoarseGrid.DifferentiationPaths{tp}(:,1)),...
%     expdata.CoarseGrid.y( expdata.CoarseGrid.DifferentiationPaths{tp}(:,2)),...
%     'ok','MarkerSize',4,'MarkerFaceColor','w');


pos = get(gca,'Position');
pos(2) = pos(2)+0.1;
pos(4) = pos(4)-0.1;
set(gca,'Position',pos);
if(tp==1), ylabel(expdata.CoarseGrid.ychanName);end
if(tp==1),
    %                 set(gca,'Units','Inches');
    c=colorbar();
    c.Location = 'EastOutside';
    %                 c.Position(1) = c.Position(1)+0.5;
    c.Ticks = [];
    set(gca,'Position',pos);
end
xlabel(expdata.CoarseGrid.xchanName);
     
% Get PCA of raw data in log space and overlay first component on heatmap

eventsVals = [];
for ss=1:height(expdata.tabSamples)
    tempvalx = expdata.Data{ss}.(expdata.CoarseGrid.xchan);
    tempvaly = expdata.Data{ss}.(expdata.CoarseGrid.ychan);
    tempval = [tempvalx';tempvaly']';
    eventsVals = [eventsVals;  tempval];
end

eventsVals = log10(eventsVals);
eventsVals(isinf(eventsVals)) = NaN;
cv = nancov(eventsVals);
[vec,val] = eig(cv);

[~,mx] = max(diag(val));
largestvec = vec(:,mx);
densityout.angl = atan(largestvec(2)/largestvec(1));


densityout.mean = nanmean(eventsVals);
densityout.median = nanmedian(eventsVals);

if(toPlotPCA)
    len = 3;
    plot(densityout.mean(1)+[0,cos(densityout.angl)*len],densityout.mean(2)+[0,sin(densityout.angl)*len],'--k','LineWidth',1.5);
    plot(densityout.mean(1)-[0,cos(densityout.angl)*len],densityout.mean(2)-[0,sin(densityout.angl)*len],'--k','LineWidth',1.5);    
end

% plot(densityout.mean(1),densityout.mean(2),'*w','MarkerSize',8);

