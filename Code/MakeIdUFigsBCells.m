% 1) Imports txt data
% 2) ScatterSlices flexibly
% 3) Fits pairwise 
% 4) Analyzes with cell-magnet to get force term
% 5) Correlate force to differentiation velocity

Actions = categorical({'ImportTxtToExpdata','GateRelevantCells','ScatterSlice','PlotTimeSeries','PrepareForPairwise','RunFitPairwise','CalcForce','CorrForceVelocity'});
% Actions = categorical({'ScatterSlice','GateRelevantCells','PrepareForPairwise','RunFitPairwise','CalcForce','CorrForceVelocity'});


%% Import
% Main WT C57BL6 timeseries acquired Dec 2016 in NCI
dirname = '../Data/Raw/';
savedir = '../Data/;
savefile = [dirname 'RawBCells'];
Cols = 'brg';
Markers = 'os^'; 


if(any(Actions=='ImportTxtToExpdata'))
    warning('off', 'MATLAB:table:ModifiedVarnames');
   expdata = ImportCytofTxtToExpdata( [savefile]);
   expdata.Time = unique(expdata.tabSamples.Time(expdata.tabSamples.Time>=0));
   expdata.nTime = length(expdata.Time);
   expdata.Replicate = unique(expdata.tabSamples.Replicate(expdata.tabSamples.Time>=0));
   expdata.nReplicate = length(expdata.Replicate);
   expdata.Type = categorical(unique(expdata.tabSamples.Type(expdata.tabSamples.Time>=0)));
   expdata.nType = length(expdata.Type);
   expdata.savefile = savefile;
   save([expdata.savefile],'expdata');
else
    load(savefile);
end
unsliced_expdata = expdata;

%% GateRelevantCells
if(any(Actions=='GateRelevantCells'))
    % Gate is the intersect of a series of ranges

%     CoarseGrid = struct;
%     CoarseGrid.xchan = 'Sm149Di_149Sm_CD19';
%     CoarseGrid.xchanName = 'CD19';
    CoarseGrid.xchan = 'Eu151Di_151Eu_IgM';
    CoarseGrid.xchanName = 'IgM';
    
    
%     CoarseGrid.ychan = 'Gd160Di_160Gd_B220';
%     CoarseGrid.ychanName = 'B220';
%      CoarseGrid.ychan = 'Nd146Di_146Nd_CD43';
%      CoarseGrid.ychanName = 'CD43';
  CoarseGrid.ychan = 'Nd150Di_150Nd_IgD';
     CoarseGrid.ychanName = 'IgD';
    CoarseGrid.EdUChannel = 'I127Di_127I_IdU';
    %     expdata.EdUChannel = 'x_APC_A__EdU';
    
 
    % MHCII direction
    warning('Overwriting x slice size');
    CoarseGrid.SliceSizex = 0.2;
    CoarseGrid.xlim = [0 3];
    CoarseGrid.lenx = 1+ceil((max(CoarseGrid.xlim)-min(CoarseGrid.xlim))/CoarseGrid.SliceSizex);
    CoarseGrid.x = linspace(CoarseGrid.xlim(1),CoarseGrid.xlim(2),CoarseGrid.lenx);
    
    % Ly6C Direction
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

    Gates.Ranges = [10.^[1 3];10.^[1 3]]; % Not bad fit
    
    
    
%     Gates.Ranges = [10.^[1 4];10.^[2 5]]; % Not bad fit
    %     Gates.Ranges = 10.^[2.5 2.7]; % 1d trajectory similar to bio paper
    
    
% unsliced_expdata = GateRelevantCellsOld(unsliced_expdata,Gates);
% unsliced_expdata = OptManifoldGateRelevantCells(unsliced_expdata,Gates);

 figure;
 set(gcf,'Color','w');
 set(gcf,'InvertHardcopy','off');
 set(gcf,'Units','Inches');
 set(gcf, 'PaperSize', [7 2]);
 set(gcf,'Position',[3,4,7,2]);
 set(gcf,'PaperUnits','Inches');
 set(gcf, 'PaperPositionMode', 'manual');
 set(gcf,'PaperPosition',get(gcf,'Position'));
 
 unsliced_expdata = GateRelevantCells(unsliced_expdata,Gates);
    
 figname = 'BCellGating';
 %     figname = 'Fig2';
 set(gcf,'Name',figname);
%  print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
 print(gcf,'-dpng',['Fig' filesep figname],'-r300');

end

%% Scatterslice
if(any(Actions=='ScatterSlice'))
    expdata = unsliced_expdata;    
    
    % 20161209: NIH No bead depletion, with CXCR4, 59 age-matched female C57BL6
           
    % Now slice
    expdata = ScatterSlice(expdata,CoarseGrid);  
%     save(expdata.savefile,'expdata');
end


%% Plot time-series with peak time colored
if(any(Actions=='PlotTimeSeries'))
    
    vert = repmat((1:expdata.CoarseGrid.leny),expdata.CoarseGrid.lenx,1);
    vert = vert(:);
    hor = repmat((1:expdata.CoarseGrid.lenx)',expdata.CoarseGrid.leny,1);
    hor = hor(:);
    allpairs = [hor vert];
    Paths = {allpairs,allpairs,allpairs}';
    tempdata = PrepareForPairwise(expdata, Paths);
    figure;
    set(gcf,'Color','w');
    set(gcf,'InvertHardcopy','off');
    set(gcf,'Units','Inches');
    set(gcf, 'PaperSize', [7 2]);
    set(gcf,'Position',[3,4,7,2]);
    set(gcf,'PaperUnits','Inches');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperPosition',get(gcf,'Position'));  
    FullData=PlotTimeSeries(tempdata,10);

    
    disp('Recorded peak times');
    
    
    figname = 'BCellTimeseries';
%     figname = 'Fig1';
    set(gcf,'Name',figname);
    print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
    print(gcf,'-dpng',['Fig' filesep figname],'-r600');
%      print(gcf,'-depsc2',['Fig' filesep figname]);
    

  
end

%% Peak times boxplot


figure;
set(gcf,'Color','w');
set(gcf,'InvertHardcopy','off');
set(gcf,'Units','Inches');
set(gcf, 'PaperSize', [7 1.5]);
set(gcf,'Position',[3,4,7,1.5]);
set(gcf,'PaperUnits','Inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperPosition',get(gcf,'Position'));
doneLegend = zeros(FullData.nTypes,1);

subplot(1,2,1);


for tp=1:FullData.nTypes
%    subplot(1,FullData.nTypes,tp);
   set(gca,'FontSize',8);
   
   HistInType = zeros(FullData.CoarseGrid.leny,FullData.ntime);
   EventsInTime = cell(FullData.ntime,1);
   EventsInLy6C = cell(FullData.CoarseGrid.leny,1);
   
   
   MarkerSizes = FullData.CoarseGrid.ScatterSliceCellCounts{tp};
   MarkerSizes(MarkerSizes==0) = NaN;
   MarkerSizes = log2(MarkerSizes/max(MarkerSizes(:)));
   MarkerSizes(MarkerSizes<-7) = 0.25;
   MarkerSizes(MarkerSizes<-3) = 3;
   MarkerSizes(MarkerSizes<=0) = 6;
   
   
   for ss=1:height(FullData.tabPopulations)
      if(FullData.tabPopulations.Type(ss) ~= FullData.Types(tp)),
          continue;
      end
      yId = FullData.tabPopulations.yId(ss);
      xId = FullData.tabPopulations.xId(ss);
      yval = FullData.CoarseGrid.y(yId)+ (rand()-0.5)*0.1;
      XdUPeakTime = FullData.tabPopulations.XdUPeakTime(ss)+(tp-2)*4;
      if(isnan(XdUPeakTime)), continue; end
      tInd = find(FullData.time==FullData.tabPopulations.XdUPeakTime(ss),1);      
      HistInType(yId,tInd) = HistInType(yId,tInd) + 1;
      
      % Events for box plot
      if(isempty(EventsInTime{tInd})), EventsInTime{tInd} = []; end
      if(isempty(EventsInLy6C{yId})), EventsInLy6C{yId} = []; end
      EventsInTime{tInd} = [EventsInTime{tInd} 10^FullData.CoarseGrid.y(yId)];
      EventsInLy6C{yId} = [EventsInLy6C{yId} FullData.time(tInd)];
      
      XdUPeakVal = FullData.tabPopulations.XdUPeakVal(ss);
      
%       if(~isnan(MarkerSizes(yId,xId)))
          pl = plot(XdUPeakTime,10^yval,'.','Marker',Markers(tp),'MarkerSize',MarkerSizes(yId,xId),...
              'Color',Cols(tp),'MarkerFaceColor',Cols(tp),'DisplayName',char(FullData.Types(tp)));
          pl.Color(4) = 0.2;
%       end
      
      if(doneLegend(tp)==1)
          set(get(get(pl,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
      else
          doneLegend(tp) = 1;
      end
      hold on
   end

%    % Add Ly6C box plots
%    for yy=1:data.CoarseGrid.leny
%        if(isempty(EventsInLy6C{yy})), continue; end
%        bx = boxplot(EventsInLy6C{yy},'Orientation','horizontal','Positions',data.CoarseGrid.y(yy)+(tp-2)*0.4,'Color',Cols(tp),'Widths',3.5,'Notch','on','BoxStyle','filled');
%        h = findobj(bx,'tag','Box');
%        h.LineWidth = 1;
%        h=findobj(bx,'tag','Outliers');
%        h.Color = Cols(tp);
%        h.MarkerFaceColor = Cols(tp);
%        %        set(bx,'linew',1.5);
%    end

   
   % Add time box plots
   for tt=1:FullData.ntime
       if(isempty(EventsInTime{tt})), continue; end
       bx = boxplot(EventsInTime{tt},'Positions',FullData.time(tt)+(tp-2)*3,'Color',Cols(tp),'Widths',3,'PlotStyle','compact');%,'Notch','on');%,'BoxStyle','filled');
       h = findobj(bx,'tag','Box');
       h.LineWidth = 1.5;
       h=findobj(bx,'tag','Outliers'); 
       h.Color = Cols(tp);
       h.MarkerFaceColor = Cols(tp);
%        set(bx,'linew',1.5);
   end
   
   set(gca,'YScale','Log');
%    ylim([FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]); 
   ylim(10.^[FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]); 
   

end

% xtickInds = find((FullData.time<=100).*(FullData.time>0));
set(gca,'XTick',FullData.time(2:9));
set(gca,'XTickLabelMode','auto');


% xlabel('Peak Time (hr)');
xlim([0 100]);
ylabel(FullData.CoarseGrid.ychanName);
title('EdU Peak Time (hr)');

% l = legend('show');

  
% figname = 'PeakTimeBoxPlot';
% %     figname = 'Fig1';
% set(gcf,'Name',figname);
% print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
% print(gcf,'-dpng',['Fig' filesep figname],'-r600');
%     
% Peak vals boxplot

subplot(1,2,2);

% figure;
% set(gcf,'Color','w');
% set(gcf,'InvertHardcopy','off');
% set(gcf,'Units','Inches');
% set(gcf, 'PaperSize', [3.42 1.5]);
% set(gcf,'Position',[3,4,3.42 1.5]);
% set(gcf,'PaperUnits','Inches');
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperPosition',get(gcf,'Position'));



doneLegend = zeros(FullData.nTypes,1);

MarkerSizes = FullData.CoarseGrid.ScatterSliceCellCounts{tp};
MarkerSizes(MarkerSizes==0) = NaN;
MarkerSizes = log2(MarkerSizes/max(MarkerSizes(:)));
MarkerSizes(MarkerSizes<-7) = 0.25;
MarkerSizes(MarkerSizes<-3) = 3;
MarkerSizes(MarkerSizes<=0) = 6;

for tp=1:FullData.nTypes
%    subplot(1,FullData.nTypes,tp);
   set(gca,'FontSize',8);
   
   HistInType = zeros(FullData.CoarseGrid.leny,FullData.ntime);
   EventsInTime = cell(FullData.ntime,1);
   EventsInLy6C = cell(FullData.CoarseGrid.leny,1);
   
   
   for ss=1:height(FullData.tabPopulations)
      if(FullData.tabPopulations.Type(ss) ~= FullData.Types(tp)),
          continue;
      end
      
      yId = FullData.tabPopulations.yId(ss);
      xId = FullData.tabPopulations.xId(ss);
      yval = FullData.CoarseGrid.y(yId)+ (tp-2)*0.25;
      XdUPeakVal = FullData.tabPopulations.XdUPeakVal(ss);
%       XdUPeakTime = FullData.tabPopulations.XdUPeakTime(ss)+(tp-2)*3;
      if(isnan(XdUPeakVal)), continue; end
%       tInd = find(FullData.time==FullData.tabPopulations.XdUPeakTime(ss),1);      
%       HistInType(yId,tInd) = HistInType(yId,tInd) + 1;
      
      % Events for box plot
      if(isempty(EventsInTime{tInd})), EventsInTime{tInd} = []; end
      if(isempty(EventsInLy6C{yId})), EventsInLy6C{yId} = []; end
%       EventsInTime{tInd} = [EventsInTime{tInd} 10^FullData.CoarseGrid.y(yId)];
      EventsInLy6C{yId} = [EventsInLy6C{yId} XdUPeakVal];
      
      
      
      if(~isnan(MarkerSizes(yId,xId)))
          pl = plot(XdUPeakVal,10^yval,'.','Marker',Markers(tp),'MarkerSize',MarkerSizes(yId,xId),...
          'Color',Cols(tp),'MarkerFaceColor',Cols(tp),'DisplayName',char(FullData.Types(tp)));
      end
      
      if(doneLegend(tp)==1)
          set(get(get(pl,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
      else
          doneLegend(tp) = 1;
      end
      hold on
   end

%    % Add Ly6C box plots (horizontal)
% for yy=1:data.CoarseGrid.leny
%     if(isempty(EventsInLy6C{yy})), continue; end
%     bx = boxplot(EventsInLy6C{yy},'Orientation','horizontal','Positions',10^(data.CoarseGrid.y(yy)+(tp-2)*0.25),'Color',Cols(tp),'Widths',3.5);%,'PlotStyle','Compact');%'Notch','on','BoxStyle','filled');
%     h = findobj(bx,'tag','Box');
%     h.LineWidth = 1;
%     h=findobj(bx,'tag','Outliers');
%     h.Color = Cols(tp);
%     h.MarkerFaceColor = Cols(tp);
%     %        set(bx,'linew',1.5);
% end
  
   % Add time box plots
%    for tt=1:data.ntime
%        if(isempty(EventsInTime{tt})), continue; end
%        bx = boxplot(EventsInTime{tt},'Positions',data.time(tt)+(tp-2)*3,'Color',Cols(tp),'Widths',3.5,'Notch','on','BoxStyle','filled');
%        h = findobj(bx,'tag','Box');
%        h.LineWidth = 1;
%        h=findobj(bx,'tag','Outliers'); 
%        h.Color = Cols(tp);
%        h.MarkerFaceColor = Cols(tp);
% %        set(bx,'linew',1.5);
%    end
   

end

set(gca,'YScale','Log');
ylabel(FullData.CoarseGrid.ychanName);
set(gca,'ytickmode','auto')
set(gca,'yticklabelmode','auto')
ylim(10.^[FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]);

 
% set(gca,'XScale','Log');
%    ylim([FullData.CoarseGrid.y(1) FullData.CoarseGrid.y(end)]);
% set(gca,'XTick',10.^[0:4]);
% set(gca,'XTickLabelMode','auto');
% xlim(10.^[0,4]);
xlim([0 400]);
set(gca,'XTick',[0:100:400]);
xlabel('EdU MFI');

if(FullData.nTypes>1)
    l = legend('show');
    l.Location = 'NorthWest';
    l.Units = 'Inches';
    l.Interpreter = 'None';
    pos = l.Position;
    pos(1) = pos(1) - 0.05;
    pos(2) = pos(2) + 0.1;
    l.Position = pos;
end

title('EdU Peak MFI');


figname = 'BCellPeakValBoxPlot';
%     figname = 'Fig1';
set(gcf,'Name',figname);
print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
print(gcf,'-dpng',['Fig' filesep figname],'-r600');
       

%% PrepareForPairwise - collect populations, average replicates, construct time-series
if(any(Actions=='PrepareForPairwise'))
    % Choose only MHCII- cells (xId==1), differentiate along Ly6C direction (yId) high to low 
%     DifferentiationPath = [ones(expdata.CoarseGrid.leny,1),(expdata.CoarseGrid.leny:-1:1)']; % xId, yId
%     DifferentiationPath =[4*ones(expdata.CoarseGrid.leny,1),(expdata.CoarseGrid.leny:-1:1)']; % xId, yId

    expdata.CoarseGrid.DifferentiationPaths = cell(expdata.nType,1);
    
    % Hardcode differentiation path: messy
    % BM:
%     expdata.CoarseGrid.DifferentiationPaths{1} = [4,13;4,12;4,11;4,10;4,9;4,8;4,7;4,6;4,5;4,4;4,3];
%     % Blood:
%     %      expdata.CoarseGrid.DifferentiationPaths{2} = [5,14;5,13;5,12;5,11;5,10;5,9;5,8;5,7;4,6;4,5;4,4;4,3];
%     expdata.CoarseGrid.DifferentiationPaths{2} = [5,13;5,12;5,11;5,10;5,9;5,8;5,7;4,6;4,5;4,4;4,3];
%     % Spleen:
%     %     expdata.CoarseGrid.DifferentiationPaths{3} = [5,14;5,13;5,12;5,11;5,10;5,9;5,8;5,7;5,6;5,5;5,4;5,3];
%     expdata.CoarseGrid.DifferentiationPaths{3} = [5,13;5,12;5,11;5,10;5,9;5,8;5,7;5,6;5,5;5,4;5,3];
%     
    
    

    % Plot each tissue in scatterslice form and show differentiation path
    % Then prepare data for pairwise
    figure;
%     miny = 2.4;
%     maxy = 4.4;
    miny = 0;
    maxy = 3;
    
    for tp=1:length(expdata.Type)

%         % Find max for each y
%         [m,mind] = max(expdata.CoarseGrid.ScatterSliceCellCounts{tp},[],2);
%         fn = find((expdata.CoarseGrid.y<=maxy).*(expdata.CoarseGrid.y>=miny));
%         fn = fn(end:-1:1);
%         expdata.CoarseGrid.DifferentiationPaths{tp} = zeros(length(fn),2);
%         expdata.CoarseGrid.DifferentiationPaths{tp}(:,1) = mind(fn);
%         inds = 1:expdata.CoarseGrid.leny;
%         expdata.CoarseGrid.DifferentiationPaths{tp}(:,2) = inds(fn);

        % Hardcode path
%         temppath = [7 7; 8 7; 9 7; 10 7; 11 7 ; 12 7; 13 7; 14 7; ...
%             14 8; 14 9; 14 10; 14 11; 13 11; 12 11 ; 11 11; 10 11; ...
%             9 11; 8 11; 7 11; 7 12; 7 13];
          temppath = [8 7; 8 8; 8 9; 8 10; 8 11; 8 12; 8 13];
%         temppath = [18 18;17 17;16 16;15 15; 14 14;13 13;12 12;11 11; 10 10; 9 9; 8 8; 7 7; 6 6; 5 5; 4 4];
        expdata.CoarseGrid.DifferentiationPaths{tp} = temppath;

        to_plot = true;
        if(to_plot)
            subplot(1,expdata.nType,tp);
            imagesc(expdata.CoarseGrid.x, expdata.CoarseGrid.y,...
                log10(expdata.CoarseGrid.ScatterSliceCellCounts{tp}));
            hold on
            set(gca,'YDir','Normal');
            colormap(bone(256));
            set(gca,'FontSize',7);
            t=title(char(expdata.Type(tp)));
            t.Interpreter = 'None';
            set(gcf,'Name',['ScatterSlice Density ' char(expdata.Type(tp))]);
            
%             ytick = expdata.CoarseGrid.y(1:2:end);
            ytick = [1:0.5:3];
            set(gca,'YTick',ytick);
            yticklab = cell(length(ytick),1);
            for yt=1:length(yticklab)
                yticklab{yt} = ['10^{' num2str(ytick(yt)) '}'];
            end
            set(gca,'YTickLabel',yticklab);
               
            xtick = [1:0.5:3];
            set(gca,'XTick',xtick);
            xticktab = cell(length(xtick),1);
            for xt=1:length(xticktab)
                xticktab{xt} = ['10^{' num2str(xtick(xt)) '}'];
            end
            set(gca,'XTickLabel',xticktab);
            xlim([1 3]);
%             set(gca,'XTick',expdata.CoarseGrid.x(1:4:end));
%                     set(gca,'XTickLabel',cellstr(['10^{' num2str(expdata.CoarseGrid.x(1:4:end)) '}']));
%             
            plot(expdata.CoarseGrid.x( expdata.CoarseGrid.DifferentiationPaths{tp}(:,1)),...
                expdata.CoarseGrid.y( expdata.CoarseGrid.DifferentiationPaths{tp}(:,2)),...
                'ok','MarkerSize',2,'MarkerFaceColor','w');
            
              
            pos = get(gca,'Position');
            pos(2) = pos(2)+0.1;
            pos(4) = pos(4)-0.1;
            set(gca,'Position',pos);
            if(tp==1), ylabel(CoarseGrid.ychanName);end
            if(tp==1),
%                 set(gca,'Units','Inches');
                c=colorbar();
                c.Location = 'EastOutside';
%                 c.Position(1) = c.Position(1)+0.5;
               set(gca,'Position',pos);
            end
             xlabel(CoarseGrid.xchanName);
            set(gca,'XTickLabelRotation',45);     
        end
        
    end
    
    set(gcf,'Color','w');
    set(gcf,'InvertHardcopy','off');
    set(gcf,'Units','Inches');
    set(gcf, 'PaperSize', [3.42 2]);
    set(gcf,'Position',[3,4,3.42,2]);
    set(gcf,'PaperUnits','Inches');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperPosition',get(gcf,'Position'));
    
    figname = 'BCellDifferentiationPath';
%     figname = 'Fig1';
    set(gcf,'Name',figname);
%     print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
    print(gcf,'-dpng',['Fig' filesep figname],'-r300');
    
    
    % Get differentiation paths for max density in x direction, down y direction
%     DifferentiationPath = MaxDensityXDifferentiationPath(expdata);
    

    data = PrepareForPairwise(expdata, expdata.CoarseGrid.DifferentiationPaths);   
    data.filename = [expdata.dirname 'TimeseriesData'];
    save( data.filename,'data');

end

% %% Plot CXCR4+Ly6C+ to BM and blood
% 
% tabMFI = readtable('CXCr4+Ly6C+EdUMFI.csv');
% tabJoined = join(tabMFI, expdata.tabSamples);
% 
% CXCR4val = zeros(expdata.nTime,1);
% CXCR4val_err = zeros(expdata.nTime,1);
% for tt=1:expdata.nTime
%     fn = find(tabJoined.Time==expdata.Time(tt));
%     if(isempty(fn)), continue; end
%     CXCR4val(tt) = mean(tabJoined.EdUMFI(fn));
%     CXCR4val_err(tt) = std(tabJoined.EdUMFI(fn))/sqrt(length(fn));
% end
% CXCR4val = CXCR4val - min(CXCR4val);
% 
% figure;
% set(gcf,'Color','w');
% set(gcf,'InvertHardcopy','off');
% set(gcf,'Units','Inches');
% set(gcf, 'PaperSize', [7 3]);
% set(gcf,'Position',[3,4,7,3]);
% set(gcf,'PaperUnits','Inches');
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperPosition',get(gcf,'Position'));
% 
% deltaType = data.nPopulations/data.nTypes;
% for tp=1:data.nTypes
%     for dd=1:6
%         subplot(3,6,dd+(tp-1)*6);
%         set(gca,'FontSize',8);
%         hold on        
%         destpop = (tp-1)*deltaType+dd;
%         ftime = find((data.time>0).*(data.time<140));
%         errorbar(data.time(ftime),CXCR4val(ftime),CXCR4val_err(ftime),'o-b','MarkerSize',3,'MarkerFaceColor','b','DisplayName','Source: CXCR4+Ly6C+');
%         
%         
% %         fitout = fitPairParseArgs(data.time(ftime),CXCR4val(ftime),data.MFI_all(ftime,destpop),data.MFI_all_err(ftime,destpop),'BruteForce',true);
%         fitout = fitPairParseArgs(data.time(ftime),CXCR4val(ftime),data.MFI_all(ftime,destpop),data.MFI_all_err(ftime,destpop),'BruteForce',false);
%         errorbar(data.time(ftime),data.MFI_all(ftime,destpop),data.MFI_all_err(ftime,destpop),'.r','Marker',Markers(tp),'MarkerSize',3,'MarkerFaceColor','r');
%         plot(data.time(ftime),fitout.ypred,'-r','LineWidth',1);
% %         set(gca,'YScale','Log');
%         dest_fit_max = fitout.ypred+fitout.delta;
%         dest_fit_min = fitout.ypred-fitout.delta;
%         
%         h = area(data.time(ftime),[dest_fit_min'; dest_fit_max'-dest_fit_min']','LineStyle','None');
%         h(1).FaceAlpha = 0;
%         h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
%         h(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
%         h(2).FaceAlpha = 0.25;
%         h(2).FaceColor = [0 0 0];
% 
% 
% 
%         ylim([0 6000]);
%         xlim([0 max(data.time(ftime))]);
%         
%         % Verify index
%         tempind = data.tabPopulations.Index(destpop);
%         if(tempind~=destpop)
%             error('Struct data mismatch !');
%         end
%         log10Ly6C = data.CoarseGrid.y(data.tabPopulations.yId(destpop));
%         if(dd==1)
%             t=title([char(data.Types(tp)) ' Ly6C 10^{' num2str(data.CoarseGrid.y(data.tabPopulations.yId(destpop))) '}']);
%         else
%             t=title(['Ly6C 10^{' num2str(data.CoarseGrid.y(data.tabPopulations.yId(destpop))) '}']);
%         end
% %         t.Units = 'Inches';
%         t.FontSize = 8;
% %         pos = t.Position;
% %         pos(2) = pos(2) - 0.2;
% %         pos(1) = pos(1) + 0.15;
% %         t.Position = pos;
%         
%         if(tp==data.nTypes)
%             xlabel(['Time (hr)']);
%         else
%             xlabel('');
%         end
%         if(dd==1)
%             ylabel('EdU MFI');
%             set(gca,'YTick',[0:2500:5000]);
%         else
%             ylabel('');
%             set(gca,'YTick',[]);
%         end
% %         grid;
%         set(gca,'XTick',[0 data.time(ftime(end))]);
%         
%         
%     end
% end
% 
% figname = 'CXCR4SrcDestFit';
% %     figname = 'Fig3';
% set(gcf,'Name',figname);
% print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
% print(gcf,'-dpng',['Fig' filesep figname],'-r600');
% %      print(gcf,'-depsc2',['Fig' filesep figname]);
% 
% 

%% FitPairwise - Run pairwise transfer function 
%  fit pairs between nearest neighbors to extract velocities
if(any(Actions=='RunFitPairwise'))
    LinearPathIndeces = zeros(length(expdata.CoarseGrid.DifferentiationPaths{1}),data.nTypes);
    for tp=1:data.nTypes
        f = find(data.tabPopulations.Type==data.Types(tp),1);
        LinearPathIndeces(:,tp) = f:(f+length(expdata.CoarseGrid.DifferentiationPaths{tp})-1);
    end
    LinearPathIndeces(:,1) = LinearPathIndeces(end:-1:1,1);
    fitouts = RunFitPairwise(data,LinearPathIndeces,data.CoarseGrid.phenmap);    
    set(gca,'YScale','Linear');
    ylim([0 300]);
    
    set(gcf,'Color','w');
    set(gcf,'InvertHardcopy','off');
    set(gcf,'Units','Inches');
    set(gcf,'PaperUnits','Inches');
    set(gcf, 'PaperSize', [7 2]);
    set(gcf,'Position',[0,0,7,2]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperPosition',get(gcf,'Position'));
    
    figname = 'BCellTimeseriesFit';
%     figname = 'Fig1';
    set(gcf,'Name',figname);
%     print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
    tic;
    print(gcf,'-dpng',['Fig' filesep figname],'-r300');
    toc;
    
    
    
end

%% Now correlate the flux in with the flux out

Markers = 'os^';

RealKappas = fitouts.kappas;
RealKappasErr = fitouts.kappas_err;
RealDs = fitouts.ds;
RealDsErr = fitouts.ds_err;
%
figure;
figname = 'BCellTotalFlux';
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

%% Quit now

disp('Finished basic plots');
gggg;

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


