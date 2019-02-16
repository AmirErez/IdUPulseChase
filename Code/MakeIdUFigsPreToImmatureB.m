% Make Neutrophil plots for IdU paper

%% Import
dirname = '../Data/';
savefile = [dirname 'RawPreBToImmatureB'];
Cols = 'brg';
Markers = 'os^'; 

warning('off', 'MATLAB:table:ModifiedVarnames');
expdata = ImportCytofTxtToExpdata( [savefile]);
expdata.Time = unique(expdata.tabSamples.Time(expdata.tabSamples.Time>=0));
expdata.nTime = length(expdata.Time);
expdata.Replicate = unique(expdata.tabSamples.Replicate(expdata.tabSamples.Time>=0));
expdata.nReplicate = length(expdata.Replicate);
expdata.Type = categorical(unique(expdata.tabSamples.Type(expdata.tabSamples.Time>=0)));
expdata.nType = length(expdata.Type);
expdata.savefile = savefile;
save(expdata.savefile,'expdata');

%     load(savefile);
unsliced_expdata = expdata;

%% Stability analysis
tabStability = readtable('../Data/RawPreBToImmatureB/StabilityPreToImmature.csv');

tabStability = join(expdata.tabSamples,tabStability);

figname = 'PreBToImmatureStability';
newfigure(figname,3.42,0.8);

for cc=1:4
   subplot(1,4,cc);
   hold on
   tableColumn = cc+width(tabStability)-4;
   plot(tabStability.Time,tabStability{:,tableColumn},'.');
   t=title(tabStability.Properties.VariableNames{tableColumn});
   t.FontSize = 9;
   xlim([min(tabStability.Time),    max(tabStability.Time)]);
   if(cc==1)
       ylim([0 100]);
%        xlabel('Time (hr)');
   else
       ylim([0.1 1000]);
       set(gca,'yscale','log');
       set(gca,'YTick',10.^[-1 1 3]);
   end
   
end

print(gcf,'-dpng',['Fig' filesep figname],'-r600');
%     print(gcf,'-dtiff',['Fig' filesep figname],'-r600');


%% GateRelevantCells

% Gate is the intersect of a series of ranges

%     CoarseGrid = struct;
CoarseGrid.xchan = 'Nd146Di_146Nd_CD43';
CoarseGrid.xchanName = 'CD43';
CoarseGrid.ychan = 'Eu151Di_151Eu_IgM';
CoarseGrid.ychanName = 'IgM';
CoarseGrid.EdUChannel = 'I127Di_127I_IdU';
    
 
% x direction
warning('Overwriting x slice size');
CoarseGrid.SliceSizex = 0.2;
CoarseGrid.xlim = [0 2];
CoarseGrid.lenx = 1+ceil((max(CoarseGrid.xlim)-min(CoarseGrid.xlim))/CoarseGrid.SliceSizex);
CoarseGrid.x = linspace(CoarseGrid.xlim(1),CoarseGrid.xlim(2),CoarseGrid.lenx);
    
% y Direction
warning('Overwriting y slice size');
CoarseGrid.SliceSizey = 0.2;
CoarseGrid.ylim = [0 3];
CoarseGrid.leny = 1+ceil((max(CoarseGrid.ylim)-min(CoarseGrid.ylim))/CoarseGrid.SliceSizey); % 4/(1/8) + 1 = (max-min)/dy + 1
CoarseGrid.y = linspace(CoarseGrid.ylim(1),CoarseGrid.ylim(2),CoarseGrid.leny)';
%     CoarseGrid.outdir = 'CoarseGrid16Out';
CoarseGrid.phenmap = colormap(parula(CoarseGrid.leny+2));
CoarseGrid.phenmap = CoarseGrid.phenmap(1:(end-2),:);
unsliced_expdata.CoarseGrid = CoarseGrid;
    
Gates = struct;
Gates.Channels = {unsliced_expdata.CoarseGrid.xchan,unsliced_expdata.CoarseGrid.ychan};

Gates.Ranges = [10.^[0 1.7];10.^[0 3]]; % Not bad fit
        

figname = 'PreBToImmatureGating';
newfigure(figname,7,2);
 
unsliced_expdata = GateRelevantCells(unsliced_expdata,Gates);
    
%  print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r300');

%% Scatterslice
expdata = unsliced_expdata;

% Now slice
expdata = ScatterSlice(expdata,CoarseGrid);
%     save(expdata.savefile,'expdata');


%% Plot time-series with peak time colored

vert = repmat((1:expdata.CoarseGrid.leny),expdata.CoarseGrid.lenx,1);
vert = vert(:);
hor = repmat((1:expdata.CoarseGrid.lenx)',expdata.CoarseGrid.leny,1);
hor = hor(:);
allpairs = [hor vert];
Paths = {allpairs,allpairs,allpairs}';
tempdata = PrepareForPairwise(expdata, Paths);

figname = 'PreBToImmatureTimeseries';
% newfigure(figname,7,2);
newfigure(figname,3.38,2);
FullData=PlotTimeSeries(tempdata,10);

print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r600');

%% PrepareForPairwise - collect populations, average replicates, construct time-series


figname = 'PreBToImmaturePeakTimeBoxPlot';
% figmaxima = newfigure(figname,7,1.5);
figmaxima = newfigure(figname,3.42,2);
expdata.CoarseGrid.DifferentiationPaths = cell(expdata.nType,1);    
 
% subplot(1,2,1);
    
miny = 0;
maxy = 3;    
tp=1;

%         % Find max for each y
%         [m,mind] = max(expdata.CoarseGrid.ScatterSliceCellCounts{tp},[],2);
%         fn = find((expdata.CoarseGrid.y<=maxy).*(expdata.CoarseGrid.y>=miny));
%         fn = fn(end:-1:1);
%         expdata.CoarseGrid.DifferentiationPaths{tp} = zeros(length(fn),2);
%         expdata.CoarseGrid.DifferentiationPaths{tp}(:,1) = mind(fn);
%         inds = 1:expdata.CoarseGrid.leny;
%         expdata.CoarseGrid.DifferentiationPaths{tp}(:,2) = inds(fn);

        % Hardcode path
%         temppath = [18 18;17 17;16 16;15 15; 14 14;13 13;12 12;11 11; 10 10; 9 9; 8 8; 7 7; 6 6; 5 5; 4 4];
%         expdata.CoarseGrid.DifferentiationPaths{tp} = temppath;

%         % Hardcode path for PCA about mean density
%         temppath = [20 16;19 15; 18 15;17 14;16 14;15 13; 14 13; 13 12; 12 12; 11 11; 10 11; 9 10; 8 10; 7 9; 6 9; 5 8; 4 8];

        % Hardcode path for PCA about median density
%         temppath = [20 16;19 15; 18 15;17 14;16 14;15 14; 14 13; 13 13; 12 12; 11 12; 10 11; 9 11; 8 10; 7 10; 6 10; 5 9; 4 9; 3 8];

  
ytick = [0:0.5:3];
xtick = [0:0.5:2];

toPCA = true;
densityout = PlotDensityMap(expdata,tp,ytick,xtick,false);

% Find coarse grid points along the PCA
% sPCA=[-0.65:0.01:1]; % Lengh along the PCA

sPCA=-0.9:0.01:1.6; % Lengh along the PCA

[temppath,AlongPCA] = GetpathAlongPCA(expdata,sPCA,densityout);
% cmap = colormap(parula(length(temppath)+1));

% Inverse path
temppath = temppath(end:-1:1,:);
    
expdata.CoarseGrid.DifferentiationPaths{tp} = temppath;
plot(expdata.CoarseGrid.x( expdata.CoarseGrid.DifferentiationPaths{tp}(:,1)),...
    expdata.CoarseGrid.y( expdata.CoarseGrid.DifferentiationPaths{tp}(:,2)),...
    '*k','MarkerSize',8);

% Plot width around PCA:
% plot(expdata.CoarseGrid.x(AlongPCA.shift1(:,1)),...
%     expdata.CoarseGrid.y(AlongPCA.shift1(:,2)),...
%     'ow','MarkerFaceColor','w','MarkerSize',3);
% plot(expdata.CoarseGrid.x(AlongPCA.shift2(:,1)),...
%     expdata.CoarseGrid.y(AlongPCA.shift2(:,2)),...
%     'ow','MarkerFaceColor','w','MarkerSize',3);


%     print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
%     print(gcf,'-dpng',['Fig' filesep figname],'-r300');
        
data = PrepareForPairwise(expdata, expdata.CoarseGrid.DifferentiationPaths);
data.filename = [expdata.dirname 'TimeseriesData'];
save( data.filename,'data');

print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r600');

%% Peak times plot along PCA1

if(false)
    subplot(1,2,2);
    hold on
    
    tp=1;
    set(gca,'FontSize',8);
    
    for ss=1:length(AlongPCA.sPCA)
        f = find((FullData.tabPopulations.xId==temppath(end-ss+1,1)).*...
            (FullData.tabPopulations.yId==temppath(end-ss+1,2)));
        if(length(f)~=1)
            warning('Problem')
        end
        xval = FullData.tabPopulations.XdUPeakTime(f);
        yval = AlongPCA.sPCA(ss);
        plot(xval,yval,'xb','MarkerSize',8);
        
        f = find((FullData.tabPopulations.yId==AlongPCA.shift1(ss,1)).*...
            (FullData.tabPopulations.xId==AlongPCA.shift1(ss,2)));
        if(length(f)~=1)
            warning('Problem')
        end
        xval = FullData.tabPopulations.XdUPeakTime(f);
        yval = AlongPCA.sPCA(ss);
        plot(xval,yval,'.r','MarkerSize',8);
        
        f = find((FullData.tabPopulations.yId==AlongPCA.shift2(ss,1)).*...
            (FullData.tabPopulations.xId==AlongPCA.shift2(ss,2)));
        if(length(f)~=1)
            warning('Problem')
        end
        xval = FullData.tabPopulations.XdUPeakTime(f);
        yval = AlongPCA.sPCA(ss);
        plot(xval,yval,'.r','MarkerSize',8);
        
    end
    
    ylabel('PC1');
    xlabel('EdU Peak Time (hr)');
    
    % l = legend('show');
    
    print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
    print(gcf,'-dpng',['Fig' filesep figname],'-r600');
end

%% FitPairwise - Run pairwise transfer function 

LinearPathIndeces = zeros(length(expdata.CoarseGrid.DifferentiationPaths{1}),data.nTypes);
for tp=1:data.nTypes
    f = find(data.tabPopulations.Type==data.Types(tp),1);
    LinearPathIndeces(:,tp) = f:(f+length(expdata.CoarseGrid.DifferentiationPaths{tp})-1);
end
% Reverse order
LinearPathIndeces(:,1) = LinearPathIndeces(end:-1:1,1);

figname = 'PreBToImmatureTimeseriesFit';
% newfigure(figname,7,2);
newfigure(figname,3.42,2);
cmap = colormap(parula(length(temppath)+1));
fitouts = RunFitPairwise(data,LinearPathIndeces,cmap);
set(gca,'YScale','Linear');
ylim([0 300]);

print(gcf,'-dpng',['Fig' filesep figname],'-r600');
% print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
 
% Now save rates in table
savefile = 'RatesPreToImm';
tempstruct = struct;
tempstruct.k = fitouts.kappas;
tempstruct.k_err = fitouts.kappas_err;
tempstruct.d = fitouts.ds;
tempstruct.d_err = fitouts.ds_err;
tempstruct.cell_frac = nanmean(data.CellFrac,1)';
tempstruct.cell_frac_err = nanstd(data.CellFrac_err)';
writetable(struct2table(tempstruct), [savefile '.csv']);
%% Quit now

disp('Finished basic plots');
gggg;

% %% Peak times boxplot
% 
% figure(figmaxima);
% subplot(1,2,2);
% hold on
% 
% for tp=1:FullData.nTypes
% %    subplot(1,FullData.nTypes,tp);
%    set(gca,'FontSize',8);
%    
%    HistInType = zeros(FullData.CoarseGrid.leny,FullData.ntime);
%    EventsInTime = cell(FullData.ntime,1);
%    EventsInLy6C = cell(FullData.CoarseGrid.leny,1);
%    
%    
%    MarkerSizes = FullData.CoarseGrid.ScatterSliceCellCounts{tp};
%    MarkerSizes(MarkerSizes==0) = NaN;
%    MarkerSizes = log2(MarkerSizes/max(MarkerSizes(:)));
%    MarkerSizes(MarkerSizes<-7) = 0.25;
%    MarkerSizes(MarkerSizes<-3) = 3;
%    MarkerSizes(MarkerSizes<=0) = 6;
%    
%    
%    
%    
%    for ss=1:height(FullData.tabPopulations)
%       if(FullData.tabPopulations.Type(ss) ~= FullData.Types(tp))
%           continue;
%       end
%       yId = FullData.tabPopulations.yId(ss);
%       xId = FullData.tabPopulations.xId(ss);
% %       yval = FullData.CoarseGrid.y(yId)+ (rand()-0.5)*0.1;
%       
%       % Project on main PCA axis along [0,1]
%       % Which means xmax^2+ymax^2 = 1;
%       % Normalize xId and yId accordingly,
% %       xnorm = xId / FullData.CoarseGrid.lenx / sqrt(FullData.CoarseGrid.x(end)^2+FullData.CoarseGrid.y(end)^2);
% %       ynorm = yId / FullData.CoarseGrid.leny / sqrt(FullData.CoarseGrid.x(end)^2+FullData.CoarseGrid.y(end)^2);
%       
%      
%     
%       x1=FullData.CoarseGrid.x(xId);
%       y1=FullData.CoarseGrid.y(yId);
%       if(x1<x0), continue; end
%       theta1=atan(y1/x1);
%       
%       % yval is the cos(theta)*R with theta the angle between (x0,y0),(x,y) and PC1
%       % with PC1 going through x0,y0
%       yval = cos(theta1-densityout.angl)*sqrt((x1-x0)^2+(y1-y0)^2);            
%       
%       XdUPeakTime = FullData.tabPopulations.XdUPeakTime(ss)+(tp-2)*4;
%       if(isnan(XdUPeakTime)), continue; end
%       tInd = find(FullData.time==FullData.tabPopulations.XdUPeakTime(ss),1);      
%       HistInType(yId,tInd) = HistInType(yId,tInd) + 1;
%       
%       % Events for box plot
%       if(isempty(EventsInTime{tInd})), EventsInTime{tInd} = []; end
%       if(isempty(EventsInLy6C{yId})), EventsInLy6C{yId} = []; end
%       EventsInTime{tInd} = [EventsInTime{tInd} 10^FullData.CoarseGrid.y(yId)];
%       EventsInLy6C{yId} = [EventsInLy6C{yId} FullData.time(tInd)];
%       
%       XdUPeakVal = FullData.tabPopulations.XdUPeakVal(ss);
%       
% %       if(~isnan(MarkerSizes(yId,xId)))
%           pl = plot(XdUPeakTime,10^yval,'.','Marker',Markers(tp),'MarkerSize',MarkerSizes(yId,xId),...
%               'Color',Cols(tp),'MarkerFaceColor',Cols(tp),'DisplayName',char(FullData.Types(tp)));
%           pl.Color(4) = 0.2;
% %       end
%       
%       if(doneLegend(tp)==1)
%           set(get(get(pl,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
%       else
%           doneLegend(tp) = 1;
%       end
%       hold on
%    end
%    
%    % Add time box plots
%    for tt=1:FullData.ntime
%        if(isempty(EventsInTime{tt})), continue; end
% %        bx = boxplot(EventsInTime{tt},'Positions',FullData.time(tt)+(tp-2)*3,'Color',Cols(tp),'Widths',3,'PlotStyle','compact');%,'Notch','on');%,'BoxStyle','filled');
% %        h = findobj(bx,'tag','Box');
% %        h.LineWidth = 1.5;
% %        h=findobj(bx,'tag','Outliers'); 
% %        h.Color = Cols(tp);
% %        h.MarkerFaceColor = Cols(tp);
% %        set(bx,'linew',1.5);
%    end
%    
%    set(gca,'YScale','Log');
% %    ylim([FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]); 
% %    ylim(10.^[FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]); 
%    
% 
% end
% 
% % xtickInds = find((FullData.time<=100).*(FullData.time>0));
% set(gca,'XTick',FullData.time(2:9));
% set(gca,'XTickLabelMode','auto');
% 
% 
% % xlabel('Peak Time (hr)');
% xlim([0 100]);
% ylabel(FullData.CoarseGrid.ychanName);
% title('EdU Peak Time (hr)');
% 
% % l = legend('show'); 
% 
% print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
% print(gcf,'-dpng',['Fig' filesep figname],'-r600');
%         


%% Now correlate the flux in with the flux out

Markers = 'os^';

RealKappas = fitouts.kappas;
RealKappasErr = fitouts.kappas_err;
RealDs = fitouts.ds;
RealDsErr = fitouts.ds_err;
%
figure;
figname = 'TotalFlux';
set(gcf,'Name',['Fig' filesep figname]);
set(gca,'FontSize',10);
set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf,'PaperUnits','Inches');
set(gcf, 'PaperSize', [3.42 2]);
set(gcf,'Position',[0,0,3.42,2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperPosition',get(gcf,'Position'));
    


FluxIn = NaN(size(RealDs));
FluxInErr = NaN(size(RealDs));
FluxOut = NaN(size(RealDs));
FluxOutErr = NaN(size(RealDs));
lx = NaN(size(RealDs));
ux = NaN(size(RealDs));

Ncells = nanmean(data.CellCount,1);
Ncells = reshape(Ncells,length(Ncells)/data.nTypes,data.nTypes);
Ncells(:,1) = Ncells(:,1)/nansum(Ncells(:,1))*4*10^6;
% Ncells(:,2) = Ncells(:,2)/sum(Ncells(:,2))*10^6/2;
% Ncells(:,3) = Ncells(:,3)/sum(Ncells(:,3))*10^6;


FluxOut = RealDs.*Ncells; % flux out = d*|Dest|
FluxOutErr = RealDsErr.*Ncells;
FluxIn(1:end-1,:) = RealKappas(1:(end-1),:).*Ncells(2:end,:);
FluxInErr(1:end-1,:) = RealKappasErr(1:(end-1),:).*Ncells(2:end,:);

for tp=1:data.nTypes
    allh=nan(length(RealKappas),data.nTypes); % all errorbar handles
%     erry=nanmean(abs(diff(y)));
    for k=1:length(RealKappas)
       if(isnan(FluxIn(k,tp))),continue; end
       er = errorbar(FluxIn(k,tp),FluxOut(k,tp),FluxOutErr(k,tp),Markers(tp),'Color',data.CoarseGrid.phenmap(k,:),...
            'MarkerFaceColor',data.CoarseGrid.phenmap(k,:),...
            'MarkerSize',3,'LineWidth',1,...
            'DisplayName',char(data.Types(tp)));
       hold on;
       if(k~=1),
           set(get(get(er,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
       end
        
       lx = FluxIn(k,tp)-FluxInErr(k,tp);
       ux = FluxIn(k,tp)+FluxInErr(k,tp);
       l1=line([lx ux],[FluxOut(k,tp) FluxOut(k,tp)],'Color',data.CoarseGrid.phenmap(k,:),'LineWidth',2); 
       set(get(get(l1,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
%        l2=line([lx(k) lx(k)],[y(k)-0.1*erry y(k)+0.1*erry],'Color',Cols(tp,:));
%        l3=line([ux(k) ux(k)],[y(k)-0.1*erry y(k)+0.1*erry],'Color',Cols(tp,:));
%        allh(k, 1:3)=[l1, l2, l3];
    end
%     arrayfun(@(d) set(get(get(d,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'), allh(~isnan(allh))); % exclude errorbars from legend

    set(gca,'XScale','log');
    set(gca,'YScale','log');
    hold on
    set(gca,'FontSize',10);
%     xlim([10^1 10^6]);
%     ylim([10^1 10^6]);
end

plot([10^1 10^6],[10^1 10^6],'--k','LineWidth',1,'DisplayName','in=out');
ax = gca;
ax.YTick = 10.^linspace(1,6,6);
ax.XTick = 10.^linspace(1,6,6);
ax.Units = 'Normalized';
xlabel('Cells in flux (count/hr)');
ylabel('Cells out flux (count/hr)');
l=legend('show');
set(l,'Location','northwest');
pos = get(l,'Position');
pos(1) = pos(1)-0.02;
pos(2) = pos(2)+0.04;

set(l,'Position',pos);
% Now plot t_1/2 as well
% Plot half life for differentiation in and out, which is log(2)./RealKappas or log(2)./ReadDs
ax = axes('Position',[0.63 0.3 0.3 0.25]);
for tp=1:data.nTypes
    for kk=1:(length(RealKappas)-1)
        plot(log(2)./RealKappas(kk,tp),log(2)./RealDs(kk,tp),Markers(tp),'Color',data.CoarseGrid.phenmap(kk,:),...
            'MarkerFaceColor',data.CoarseGrid.phenmap(kk,:),...
            'MarkerSize',2,...
            'DisplayName',char(data.Types(tp)));
        hold on
    end
    set(gca,'FontSize',10);
    
    maxk = max(max(log(2)./RealKappas(2:end,:)));
    maxd = max(max(log(2)./RealDs(2:end,:)));
%     xlim([0 maxk]);
%     ylim([0 maxd]);
end
grid('on');
xl=xlabel('$t_{1/2}$ in (hr)','Interpreter','Latex','FontSize',8);
set(xl,'Units','Normalized');
pos = get(xl,'Position');
pos(2) = 0.3;
set(xl,'Position',pos);
ylabel('$t_{1/2}$ out (hr)','Interpreter','Latex','FontSize',8);


print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r600');



%% Plot the linear differentiation path as "cartoon" with elipse length ~ log(abundance)

figure;
figname = 'FluxPyramid';
set(gcf,'Name',['Fig' filesep figname]);
set(gca,'FontSize',10);
set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf,'PaperUnits','Inches');
set(gcf, 'PaperSize', [7 2]);
set(gcf,'Position',[0,0,7,2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperPosition',get(gcf,'Position'));
    

for tp=1:data.nTypes
    elipseheight = 0.4;
    subplot(1,data.nTypes,tp);
    set(gca,'FontSize',10);
    hold on
    
%     norm = min(Ncells(:);
    norm = 10^(3.5);
    
    popinds = LinearPathIndeces(end:-1:1,1);
    gridvals = data.tabPopulations.yId(popinds);
    
    for yy=gridvals'
%         wdth = log10(Ncells(yy,tp));
          wdth = log10(Ncells(yy,tp)/norm*2);
          if(wdth<0), wdth = 0.1; end
          maxwidth = log10(max(Ncells(:))/norm*2);
%         wdth = Ncells(yy,tp) / min(Ncells(:));
%     maxwidth = max(Ncells(:))/min(Ncells(:));
      
     
        rectangle('Position',[-wdth/2, yy-elipseheight/2, wdth,elipseheight ],...
            'Curvature',[1,1],...
            'FaceColor',data.CoarseGrid.phenmap(yy,:));
        xa = [0 0];
        ya = [yy+0.9-elipseheight/2,yy+0.1+elipseheight/2];
        
        if(yy<max(gridvals))
%             arrowwidth = log10(RealKappas(yy,tp).*Ncells(yy-1,tp)/min(Ncells(:))*2);
%             arrowwidth = RealKappas(yy,tp).*Ncells(yy-1,tp)/min(Ncells(:))/2;
%             arrowwidth = 2;
%             arrowcol = 1-[1,1,1]*RealKappas(yy,tp).*Ncells(yy-1,tp)/max(Ncells(:));
            arrowcol = data.CoarseGrid.phenmap(yy,:);
            l = line(xa,ya,'Color',arrowcol,'LineWidth',1);
            fluxpos = find(gridvals==yy);
            fluxwidth = log10(FluxIn(fluxpos,tp)/norm*24);
            if(fluxwidth<0), fluxwidth = 0.1; end
            l = rectangle('Position',[-fluxwidth/2,ya(2),fluxwidth,ya(1)-ya(2)],'LineWidth',1);
            l = line([xa(1)-0.03 xa(2)+0.01],[ya(2)+0.15,ya(2)-0.02],'Color',arrowcol,'LineWidth',1);
            l = line([xa(1)+0.03 xa(2)-0.01],[ya(2)+0.15,ya(2)-0.02],'Color',arrowcol,'LineWidth',1);
%             [xaf,yaf] = axescoord2figurecoord(xa,ya);
%             ar = annotation('arrow',xaf,yaf,'Color',phenmap(yy,:),'LineWidth',linewidth);
        end
        
    end
    
    % Plot population size scale
     wdth = log10((10^6)/norm*2);
     crushedheight = 0.25;
     rectangle('Position',[-wdth/2, 0, wdth,crushedheight ],'LineWidth',1);%,...
 %           'Curvature',[1,1]);
     text(wdth/2 -0.1,-elipseheight*1.3,'10^6','FontSize',10);
     wdth = log10((10^5)/norm*2);
     rectangle('Position',[-wdth/2, 0, wdth,crushedheight ],'LineWidth',1);%,...
 %           'Curvature',[1,1]);
     text(wdth/2 -0.1,-elipseheight*1.3,'10^5','FontSize',10);   
     wdth = log10((10^4)/norm*2);
     rectangle('Position',[-wdth/2, 0, wdth,crushedheight ],'LineWidth',1);%,...
%            'Curvature',[1,1]);
     text(wdth/2 -0.1,-elipseheight*1.3,'10^4','FontSize',10);    

     % Labels
     rectangle('Position',[-1.3, 4,0.4,elipseheight],'LineWidth',1,'Curvature',[1,1],'FaceColor',[0 0 0.8]);
     text(-1.3,4.8,'# Cells','FontSize',10);
     rectangle('Position',[-1.3, 2,0.4,elipseheight],'LineWidth',1);
     text(-1.3,2.8,'# Cells/day','FontSize',10);
     

     
     
    xlim([-maxwidth/2*1.05 maxwidth/2*1.05]);
    set(gca,'Visible','off');
 
    an = annotation(gcf,'textbox',...
        [0.35 0.9 0.2 0.03],...
        'String',char(data.Types(tp)),...
        'LineStyle','none',...
        'FontSize',18,...
        'FitBoxToText','on');
    an.Position = [an.Position(1)-an.Position(3)/2,an.Position(2:4)];
    
end

print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r600');

%% CalcForce - Calculate force term according to cell-magnet method
if(any(Actions=='CalcForce'))
    allforces = cell(height(expdata.tabSamples),1);
    expdata.CoarseGrid.forces = zeros(length(data.CoarseGrid.DifferentiationPaths{1}),data.nTypes);
    expdata.CoarseGrid.forces_err = expdata.CoarseGrid.forces;
    
    figure;
    hold on
    xlabel('l');
    ylabel('\phi');
    set(gcf,'Name','Force by cell-magnet');
    for tp=1:expdata.nType
        for ss=1:height(expdata.tabSamples)
            % Calculate force from histogram
            if(expdata.tabSamples.Type(ss)~=expdata.Type(tp)),continue; end
            LinearPath = expdata.CoarseGrid.DifferentiationPaths{tp};
            
             % Force by savgol
             eventinds = [];
             for ll=1:length(LinearPath)
                 eventinds = [eventinds ;...
                     find((expdata.Data{ss}.xId==LinearPath(ll,1)).*...
                    (expdata.Data{ss}.yId==LinearPath(ll,2)))];
             end
             events = expdata.Data{ss}.(expdata.CoarseGrid.ychan)(eventinds);
             coarsebininds = flipud(LinearPath(:,2)); % Ly6C values are flipped to be from up to down for historgram  
             log10bins = expdata.CoarseGrid.y(coarsebininds);
             dl = log10bins(2)-log10bins(1);
             log10bins = (log10bins(1):(dl/2):log10bins(end))';
             force = CalcForceFromEvents(events, log10bins);
             force(:,1) = force(:,1)-dl/2;
            % Make histogram and calculate force by diff (noisy)
%             hst = zeros(length(LinearPath),1);
%             for ll=1:length(LinearPath)
%                 hst(ll) = nnz( (expdata.Data{ss}.xId==LinearPath(ll,1)).*...
%                     (expdata.Data{ss}.yId==LinearPath(ll,2)));
%                 
%             end
%             % Now we have the histogram, we will calculate force from it
%             l = expdata.CoarseGrid.y(LinearPath(:,2));
% %             Q = -l+log10(hst);
%             phi0 = diff(log10(hst))./diff(l)-1;            
%             force = zeros(length(l),2);
%             force(:,1) = l;
%             force(:,2) = [phi0; NaN];                       
               
            allforces{ss} = [NaN NaN; force(end:-2:2,:); NaN NaN];            
            l = plot(allforces{ss}(:,1),allforces{ss}(:,2),'-','LineWidth',2,'Color',Cols(tp),'DisplayName',char(data.Types(tp)));
            l.Color(4) = 0.1;
        end
        % Now average forces for each tissue type
        
        fn = find(expdata.tabSamples.Type==data.Types(tp));
        tempforces = zeros(length(data.CoarseGrid.DifferentiationPaths{tp}),length(fn));
        for ii=1:length(fn)
            val = allforces{fn(ii)};
            if(isempty(val)), continue; end
            tempforces(1:length(val),ii)=val(:,2);
        end
        expdata.CoarseGrid.forces(:,tp) = nanmean(tempforces,2);
        expdata.CoarseGrid.forces_err(:,tp) = nanstd(tempforces,1,2);%/sqrt(length(fn));
        plot(expdata.CoarseGrid.y(expdata.CoarseGrid.DifferentiationPaths{tp}(:,2)),...
            expdata.CoarseGrid.forces(:,tp),...
            '-','Color',Cols(tp),'Marker',Markers(tp),'MarkerFaceColor',Cols(tp),...
            'MarkerSize',12,'LineWidth',2,'DisplayName',char(data.Types(tp)));

    end
 %Now average the force per phenotype per bin to estimate errors
%     figure; hold on
%     set(gcf,'Name',['Averaged Forces']);
%     xlabel(['log_{10} ' ' ' data.CoarseGrid.ychanName]);
%     ylabel('force \phi');

   
end

%% Now correlate forces with velocities !!

% figure; hold on
% set(gcf,'Name','Velocity vs force');
% xlabel('Force');
% ylabel('Velocity');

for tp=1:data.nTypes
    figure(100+tp); 
    hold on
    set(gcf,'Name',[char(data.Types(tp)) ' Velocity vs force']);
    ylabel('Force');
    xlabel('Velocity');
    ydiffpath = expdata.CoarseGrid.DifferentiationPaths{tp}(:,2);
    for ll=1:length(ydiffpath) 
        phencolind = find(expdata.CoarseGrid.y==expdata.CoarseGrid.y(ydiffpath(ll)));
        yval = expdata.CoarseGrid.forces(ll,tp);
        yerr = expdata.CoarseGrid.forces_err(ll,tp);
        dy = expdata.CoarseGrid.y(2)-expdata.CoarseGrid.y(1);

        %         xval = fitouts.kappas(ll,tp)/dy; % Slice normalize
        xval = fitouts.kappas(ll,tp);
        xval = xval / max(fitouts.kappas(~isnan(fitouts.kappas(:,tp)),tp)); % Max normalize
        
        xvalerr = fitouts.kappas_err(ll,tp)/dy;
%         xval = fitouts.ds(ll,tp);
%         xvalerr = fitouts.ds(ll,tp);

        errorbar(xval,yval,yerr,Markers(tp),'Color',data.CoarseGrid.phenmap(phencolind,:),'MarkerSize',12,'MarkerFaceColor',data.CoarseGrid.phenmap(phencolind,:),'DisplayName',[char(expdata.CoarseGrid.ychanName) ' ' num2str(expdata.CoarseGrid.y(phencolind))]);
        hold on
        lx = xval-xvalerr;
        ux = xval+xvalerr;
        l1=line([lx ux],[yval yval],'Color',data.CoarseGrid.phenmap(phencolind,:),'LineWidth',2);
        set(get(get(l1,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
        
        
    end
%     ylim([min(data.CoarseGrid.forces(:)) max(data.CoarseGrid.forces(:))]);
%     ylim([-2.5 2.5]);
%     ylim([0 0.045]);
%     ylim([-3 1]);
    grid on
end


