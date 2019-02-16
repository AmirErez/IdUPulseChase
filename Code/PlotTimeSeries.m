% Plots time series in all tissues color coded by peak time
function [data]= PlotTimeSeries(data,maxpeaktimeind)
%% Get time of peak IdU MFI naively
tempstruct = struct;
tempstruct.XdUPeakTime = nan(data.nPopulations,1);
tempstruct.XdUPeakVal = nan(data.nPopulations,1);
tempstruct.XdUPulseWidth = nan(data.nPopulations,1);

% Define effective fit
% fitresult = cell( data.nPopulations, 1 );
% % gof = struct( 'sse', cell( data.nPopulations, 1 ), ...
% %     'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
% % ft = fittype( 'A*heaviside(x-t0)*(x-t0)^alpha*exp(-d*(x-t0))', 'independent', 'x', 'dependent', 'y' );
% % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % opts.Display = 'Off';
% % opts.Lower = [0 0 0 0];
% % opts.Upper = [Inf Inf Inf Inf];
% % opts.StartPoint = [100 1 0.1 10];
% 
% ft = fittype('poly2');
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf -Inf 0];
% opts.Upper = [0 Inf Inf];
% % opts.StartPoint = [0 0 1];
% 


fig=gcf;
fig.Units = 'Normalized';
% figs(tp) = fig;

for tp=1:data.nTypes
%  for tp=1:1

  
    % First make axes
 
           
    peakcolmap = colormap(winter(maxpeaktimeind));
    
    for pp=1:height(data.tabPopulations)
%         ax = axs{tp,data.tabPopulations.yId(pp),data.tabPopulations.xId(pp)};
%         axes(ax);

        axes('Position',[(data.tabPopulations.xId(pp)+(tp-1)*data.CoarseGrid.lenx)/(data.nTypes*data.CoarseGrid.lenx+2),...
            (data.tabPopulations.yId(pp))/(data.CoarseGrid.leny+2),...
            1/(data.CoarseGrid.lenx*data.nTypes),...
            1/(data.CoarseGrid.leny+2)])
%         ax{tp,data.tabPopulations.yId(pp),data.tabPopulations.xId(pp)} = get(gca,'Position');
        popsize = nanmean(data.CellCount(:,pp));
      
        if(data.tabPopulations.Type(pp)~=data.Types(tp) ||isnan(popsize) || popsize<25)            
            set(gca,'visible','off');
            continue;
        end
        
%         multiWaitbar('Plotting',pp/data.nPopulations);
        
        % Set peak time by naively looking for max
        [mx,mxind] = max(data.MFI_all(:,pp));
        tempstruct.XdUPeakTime(pp) = data.time(mxind);

        % Set peak by convolving with "Gaussian" 
        cnv = conv(data.MFI_all(:,pp),[0.5 1 0.5],'same');
        [mx,mxind] = max(cnv);
        tempstruct.XdUPeakTime(pp) = data.time(mxind);
        
        tempstruct.XdUPeakVal(pp) = mx;
        mxind = min(mxind,length(peakcolmap));
        PeakCol = peakcolmap(mxind,:);
        
        [x0,y0,iout,jout] = intersections(data.time,data.MFI_all(:,pp),data.time,repmat(mx/2,length(data.time),1));
        if(length(x0)==2),
            tempstruct.XdUPulseWidth(pp) = x0(2)-x0(1);
        else
            tempstruct.XdUPulseWidth(pp) = NaN;
        end
    
        % Find max by fitting effective function
        
        %     col = 'b';
        %     if(data.tabPopulations.Type(pp)=='Blood'), col = 'r'; end
        %
        f = find(data.Types==data.tabPopulations.Type(pp));
        
        %     ax = axes('Position',[(data.tabPopulations.xId(pp)-1)/data.CoarseGrid.lenx,(data.tabPopulations.yId(pp)-1)/data.CoarseGrid.leny,1/data.CoarseGrid.lenx,1/data.CoarseGrid.leny]);
        
        %         title(char(data.tabPopulations.Population(pp)));
        
        val = data.MFI_all(:,pp)-min(min(data.MFI_all));
        %     val = log10(data.MFI_all(:,pp));
        %     val_err = data.MFI_all_err(:,pp)
        
        %     val = log10(data.MFI_all(:,pp)*parameters.cell_count(pp));
        %      val =log10(data.MFI_all(:,pp)*parameters.cell_count(pp));
        %      val_err = data.MFI_all_err(:,pp)
        
        pl = plot(data.time,val,'-','Color','k','LineWidth',1);
        hold on
        
        ShiftedMFI = data.MFI_all-min(min(data.MFI_all));
        xlim([0 data.time(end)]);
%         ylim([10 max(tempval(:))]);
        
        
        h = area(data.time,[ones(size(val)) val],'LineStyle','None');
        h(1).FaceAlpha = 0;
        h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
        h(2).FaceAlpha = 1; %0.2
        h(2).FaceColor = PeakCol;
        
%         set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
        
     
        
        %     ylim([0 max([max(data.MFI_all(:,pp)),max(get(gca,'YLim'))])]);
        %     tempval = log10(data.MFI_all.*repmat(parameters.cell_count',data.ntime,1));
        
        %     ylim([min(tempval(:)) max(tempval(:))]);
        %     ylim([min(tempval(:)) max(tempval(:))])
        ylim([1 max(ShiftedMFI(:))]);
        set(gca,'XTickLabel',{});
        set(gca,'YTickLabel',{});
        
        set(gca,'XTick',[0:48:data.time(end)]);
        set(gca,'YTick',[0:0.2:4]);
        set(gca,'LineWidth',1);
        set(gca,'YScale','Log');
%         grid('on');
%         set(gca,'GridColor','w');
        set(gca,'Box','on');
        
        
        
        %     if(AlreadyPainted(data.tabPopulations.yId(pp),data.tabPopulations.xId(pp))==0)
        %         set(gca,'Color',GridColors(pp,:));
        %         AlreadyPainted(data.tabPopulations.yId(pp),data.tabPopulations.xId(pp))=1;
        %     end
        %     set(a,'LineWidth',1.5);
    end
%     axes('Position',[(data.tabPopulations.xId(pp)+(tp-1)*data.CoarseGrid.lenx)/(data.nTypes*data.CoarseGrid.lenx+2),...
%             (data.tabPopulations.yId(pp))/(data.CoarseGrid.leny+2),...
%             1/(data.CoarseGrid.lenx*data.nTypes),...
%             1/data.CoarseGrid.leny])
%         

    x1=(5+(tp-1)*data.CoarseGrid.lenx)/(data.nTypes*data.CoarseGrid.lenx+2);
    x2 = x1+1/(data.nTypes+2);
%     x2=(data.nTypes+(tp-1)*data.CoarseGrid.lenx)/(data.nTypes*data.CoarseGrid.lenx+2);
%     y1 = 1/(data.CoarseGrid.leny+2);
    y1 = 0.03;
    y2 = y1;
    annotation(gcf,'arrow',[x1 x2],...
        [y1 y2],'LineWidth',1);
    annotation(gcf,'textbox',...
        [x1-0.02 0.1 0.8/data.nTypes 0.03],...
        'String',data.CoarseGrid.xchanName,...
        'LineStyle','none',...
        'FontSize',8,...
        'FitBoxToText','off');
    
%     annotation(gcf,'textbox',...
%         [x1 0.985 0.8/data.nTypes 0.02],...
%         'String',char(data.Types(tp)),...
%         'LineStyle','none',...
%         'FontSize',10,...
%         'Interpreter','None',...
%         'FitBoxToText','off');

    % Annotate phenotype directions
    if(tp==1)
        annotation(gcf,'arrow',[0.02 0.02],...
            [0.1 0.9],'LineWidth',1);
        
        annotation(gcf,'textbox',...
            [0.02 0.8 0.2 0.03],...
            'String',data.CoarseGrid.ychanName,...
            'LineStyle','none',...
            'FontSize',8,...
            'FitBoxToText','off');
    end
    
    
    % Annotate small graphs
%     annotation(gcf,'arrow',[0.09 0.09],[0.3 0.3+1.5/data.CoarseGrid.leny],'LineWidth',1);
%     annotation(gcf,'textbox',...
%         [0.07 0.55 0.2 0.03],...
%         'String',{'MFI','IdU'},...
%         'LineStyle','none',...        
%         'FontSize',8,...
%         'FitBoxToText','off');
%     
%     annotation(gcf,'arrow',[0.09 0.09+0.5/(data.CoarseGrid.lenx*data.nTypes)],[0.3 0.3],'LineWidth',1);
%     annotation(gcf,'textbox',...
%         [0.05 0.22 0.2 0.03],...
%         'String',['Time ' num2str(data.time(end)) 'hr'],...
%         'LineStyle','none',...
%         'FontSize',8,...
%         'FitBoxToText','off');
    

    
%     % Plot master legend
%     plot(0,0,'-','Color', Cols(1,:),'LineWidth',2,'DisplayName',char(data.Types(1)));
%     hold on
%     plot(0,0,'-','Color', Cols(2,:),'LineWidth',2,'DisplayName',char(data.Types(2)));
%     plot(0,0,'-','Color', Cols(3,:),'LineWidth',2,'DisplayName',char(data.Types(3)));
    
%     l = legend('show');
%     l.Position = [0.82 0.82 0.15 0.15];
%     l.FontSize = 16;
%     l.EdgeColor = 'None';
    
end
     
% Master colorbar
if(data.nTypes>1)
    c=colorbar('Position', [1/data.nTypes-0.1  0.2 0.03 0.4]);
    annotation(gcf,'textbox',...
        [1/data.nTypes-0.1 0.87 0.05 0.03],...
        'String',{'IdU','Peak','Time', '(hr)'},...
        'LineStyle','none',...
        'FontSize',8,...
        'FitBoxToText','off');
else
    c=colorbar('Position', [0.9  0.15 0.03 0.4]);   
    annotation(gcf,'textbox',...
        [0.87 0.81 0.05 0.1],...
        'String',{'IdU','Peak','Time', '(hr)'},...
        'LineStyle','none',...
        'FontSize',8,...
        'FitBoxToText','off');    
end
c.Ticks = data.time(2:2:length(peakcolmap))/data.time(length(peakcolmap));
c.AxisLocation='out';
c.TickLabels = data.time(2:2:length(peakcolmap));

% multiWaitbar('Plotting','Close');



% Record peak time and save graphic

data.tabPopulations.XdUPeakTime = tempstruct.XdUPeakTime;
data.tabPopulations.XdUPeakVal = tempstruct.XdUPeakVal;
data.tabPopulations.XdUPulseWidth = tempstruct.XdUPulseWidth;

