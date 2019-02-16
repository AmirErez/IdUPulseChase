% Reads rates files and plots them in a single figure

filenames = {'RatesNeut.csv', 'RatesPreToImm.csv','RatesImmTransMarg.csv'};
legends = {'Neutrophils', 'Pre B to Immature B','Immature to Mature B'};
% Gate normalization: approximate fraction from total flux
GateNormalization = [0.24, 0.33*0.7, 0.33*0.25];


Markers = 'os^';
%
Cols = [0.8 0.1 0.1
        0.1 0.8 0.1
        0.1 0.1 0.8];



% figure;
% hold on
figname = 'RatesAll';
newfigure(figname,3.42,3);

% set(gcf,'Name',['Fig' filesep figname]);
% set(gca,'FontSize',10);
% set(gcf,'Color','w');
% set(gcf,'InvertHardcopy','off');
% set(gcf,'Units','Inches');
% set(gcf,'PaperUnits','Inches');
% set(gcf, 'PaperSize', [3.42 2]);
% set(gcf,'Position',[0,0,3.42,2]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperPosition',get(gcf,'Position'));

for tp=1:length(filenames)
    tabRates = readtable(filenames{tp});
    phenmap = colormap(parula(height(tabRates)+1));

    
    % Cellular flux
%     FluxIn = tabRates.k.*tabRates.cell_frac;
%     FluxInErr = tabRates.k_err.*tabRates.cell_frac;
%     FluxOut = tabRates.d(1:(height(tabRates)-1)).*tabRates.cell_frac(2:end);
%     FluxOutErr = tabRates.d_err(1:(height(tabRates)-1)).*tabRates.cell_frac(2:end);
   
 % Cellular flux as part of 1d path
%     FluxIn = tabRates.k.*tabRates.cell_frac/sum(tabRates.cell_frac);
%     FluxInErr = tabRates.k_err.*tabRates.cell_frac/sum(tabRates.cell_frac);
%     FluxOut = [tabRates.d(1:(height(tabRates)-1)).*tabRates.cell_frac(2:end)/sum(tabRates.cell_frac); NaN];
%     FluxOutErr = [tabRates.d_err(1:(height(tabRates)-1)).*tabRates.cell_frac(2:end)/sum(tabRates.cell_frac); NaN];
%    

    % Mean flux as part of 
    FluxIn = tabRates.k.*tabRates.cell_frac/nansum(tabRates.cell_frac);
    FluxInErr = NaN*FluxIn;
    FluxOut = [tabRates.d(1:(height(tabRates)-1)).*tabRates.cell_frac(2:end)/nansum(tabRates.cell_frac); NaN];
    FluxOutErr = NaN*FluxOut;
 
%     FluxIn = tabRates.cell_frac;
%     FluxOut = sqrt(tabRates.k.*tabRates.d);
    
    FluxIn = tabRates.k;
    FluxOut = tabRates.d;
    
    
%     FluxOut = sqrt(FluxIn.*FluxOut);
%     FluxIn = (0:(length(FluxIn)-1))'/(length(FluxIn)-1);
% 
%     FluxOut = sqrt(tabRates.d.*tabRates.k);
%     FluxOut = tabRates.cell_frac;
%    FluxOut = sqrt(tabRates.d.*tabRates.k);
%    FluxIn = tabRates.cell_frac;
    
    % Half life
%     FluxIn = 1./tabRates.k;
%     FluxInErr = NaN*FluxIn;
%     FluxOut = 1./tabRates.d;
%     FluxOutErr = NaN*FluxOut;
   
     plot(FluxIn,FluxOut,'*','LineWidth',1,'Color',Cols(tp,:),'DisplayName',legends{tp});

%     for k=1:length(tabRates.k)
%        if(isnan(tabRates.k(k))),continue; end
%         
%        er = errorbar(FluxIn(k),FluxOut(k),FluxOutErr(k),Markers(tp),'Color',phenmap(k,:),...
%             'MarkerFaceColor',phenmap(k,:),...
%             'MarkerSize',12,'LineWidth',1,...
%             'DisplayName',legends{tp});
%        hold on;
%        if(k~=1),
%            set(get(get(er,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
%        end
%         
%        lx = FluxIn(k)-FluxInErr(k);
%        ux = lx+2*FluxInErr(k);
%        l1=line([lx ux],[FluxOut(k) FluxOut(k)],'Color',phenmap(k,:),'LineWidth',2); 
%        set(get(get(l1,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off')
% %        l2=line([lx(k) lx(k)],[y(k)-0.1*erry y(k)+0.1*erry],'Color',Cols(tp,:));
% %        l3=line([ux(k) ux(k)],[y(k)-0.1*erry y(k)+0.1*erry],'Color',Cols(tp,:));
% %        allh(k, 1:3)=[l1, l2, l3];
%     end
% %     arrayfun(@(d) set(get(get(d,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off'), allh(~isnan(allh))); % exclude errorbars from legend

%     set(gca,'XScale','log');
%     set(gca,'YScale','log');
    hold on
%     set(gca,'FontSize',10);
%     xlim([10^-6 10^0]);
%     ylim([10^-6 10^0]);
end

% plot([10^-6 10^0],[10^-6 10^0],'--k','LineWidth',1,'DisplayName','in=out');
ax = gca;
% ax.YTick = 10.^linspace(3,6,3);
% ax.XTick = 10.^linspace(3,6,3);
ax.Units = 'Normalized';
xlabel('Cells in flux (count/hr)');
ylabel('Cells out flux (count/hr)');
l=legend('show');
set(l,'Location','northwest');
pos = get(l,'Position');
pos(1) = pos(1)-0.02;
pos(2) = pos(2)+0.04;
set(l,'Position',pos);


% % Now plot t_1/2 as well
% % Plot half life for differentiation in and out, which is log(2)./RealKappas or log(2)./ReadDs
% ax = axes('Position',[0.63 0.3 0.3 0.25]);
% for tp=1:data.nTypes
%     for kk=1:(length(RealKappas)-1)
%         plot(log(2)./RealKappas(kk,tp),log(2)./RealDs(kk,tp),Markers(tp),'Color',data.CoarseGrid.phenmap(kk,:),...
%             'MarkerFaceColor',data.CoarseGrid.phenmap(kk,:),...
%             'MarkerSize',2,...
%             'DisplayName',char(data.Types(tp)));
%         hold on
%     end
%     set(gca,'FontSize',10);
%     
%     maxk = max(max(log(2)./RealKappas(2:end,:)));
%     maxd = max(max(log(2)./RealDs(2:end,:)));
% %     xlim([0 maxk]);
% %     ylim([0 maxd]);
% end
% grid('on');
% xl=xlabel('$t_{1/2}$ in (hr)','Interpreter','Latex','FontSize',8);
% set(xl,'Units','Normalized');
% pos = get(xl,'Position');
% pos(2) = 0.3;
% set(xl,'Position',pos);
% ylabel('$t_{1/2}$ out (hr)','Interpreter','Latex','FontSize',8);


% print(gcf,'-dtiff',['Fig' filesep figname],'-r600');
% print(gcf,'-dpng',['Fig' filesep figname],'-r600');
 
