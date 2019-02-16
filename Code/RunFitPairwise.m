function fitouts = RunFitPairwise(data,LinearPathIndeces,phenmap)
% Runs linear differentiation paths as defined in data.Populations 
% Differentiation paths defined in LinearPathIndeces 
%      rows - population index in of differentiation
%      cols for different types

fitouts.kappas = NaN(length(LinearPathIndeces),data.nTypes);
fitouts.ds = NaN(length(LinearPathIndeces),data.nTypes);
fitouts.kappas_err = NaN(length(LinearPathIndeces),data.nTypes);
fitouts.SSEnormbyErr = NaN(length(LinearPathIndeces),data.nTypes);
fitouts.ds_err = NaN(length(LinearPathIndeces),data.nTypes);
fitouts.fitout = cell(length(LinearPathIndeces),data.nTypes);

% Ncells = NaN(data.CoarseGrid.leny,data.nTypes);

% figure;

fn = find(data.time>0);
data.time = data.time(fn);
data.MFI_all = data.MFI_all(fn,:);
data.MFI_all_err = data.MFI_all_err(fn,:);

for tp=1:data.nTypes
    
    subplot(1,data.nTypes,tp);
    set(gca,'FontSize',8);
    set(gcf,'Name',[char(data.Types(tp)) ' 1d propagation']);
          
    for yy=1:length(LinearPathIndeces)
        
        f=LinearPathIndeces(yy,tp);
        val = data.MFI_all(:,f) - min(data.MFI_all(:,f));
        yId = data.tabPopulations.yId(LinearPathIndeces(yy,tp));
        if(yy>1)
            pl = errorbar(data.time,val,data.MFI_all_err(:,f),'o','LineWidth',1,'Color',phenmap(yy,:),'MarkerSize',2.5,'MarkerFaceColor',phenmap(yy,:),...
                'DisplayName',[data.CoarseGrid.ychanName ' 10^{' num2str(data.CoarseGrid.y(yId),2) '}']);
            if(mod(yy,2)~=0)
                set( get( get( pl, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
        end
        hold on
%         Ncells(yy,tp) = parameters.cell_count(f);
        
        % Default error for unknown error
        data.MFI_all_err(isnan(data.MFI_all_err)) = 10;
        data.MFI_all_err(data.MFI_all_err<10) = 10;
        
        if(yy>1&&~any(isnan(prev_MFI))&&~any(isnan(val)))
%             fitout = fitPairParseArgs(data.time, prev_MFI,val, data.MFI_all_err(:,f));
            fitout = fitPairParseArgs(data.time, prev_MFI,val, data.MFI_all_err(:,f),'BruteForce',true);
            fitouts.fitout{yy-1,tp} = fitout;
            fitouts.kappas(yy-1,tp) = fitout.kappas;
            fitouts.ds(yy-1,tp) = fitout.ds;
            fitouts.kappas_err(yy-1,tp) = fitout.kappas_err;
            fitouts.ds_err(yy-1,tp) = fitout.ds_err;
            fitouts.SSEnormbyErr(yy-1,tp) = fitout.SSEnormbyErr;            
            ft = plot(data.time,fitout.ypred,'-','Color',phenmap(yy,:),'MarkerFaceColor',phenmap(yy,:),'MarkerSize',3,'LineWidth',1);
            set( get( get( ft, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            
            dest_fit_max = fitout.ypred+fitout.delta;
            dest_fit_min = fitout.ypred-fitout.delta;
            
            h = area(data.time,[dest_fit_min'; dest_fit_max'-dest_fit_min']','LineStyle','None');
            h(1).FaceAlpha = 0;
            h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
            h(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
            h(2).FaceAlpha = 0.2;
            h(2).FaceColor = phenmap(yy,:);
            
        end
        prev_MFI = val;
%         prev_MFI_err = data.MFI_all_err(:,f);
    end
    xlim([0, data.time(end)]);
    %     ylim([min(min(data.MFI_all(:,fX))),max(max(data.MFI_all(:,fX)))]);
    ylim([1 3000]);
    %     set(gca,'YScale','log');
    set(gca,'FontSize',8);
    set(gca,'YScale','log')
    xlabel('Time (hr)');
    if(tp==1), ylabel('IdU Mean Fluorescence'); end
%     t=title(char(data.Types(tp)));
%     t.Interpreter = 'None';
%     annotation(gcf,'textbox',...
%         [0.45 0.8 0.2 0.03],...
%         'String',char(data.Types(tp)),...
%         'LineStyle','none',...
%         'FontSize',22,...
%         'FitBoxToText','off');
%     if(tp==1)
%         l = legend('show');
%         l.FontSize = 6;
%         l.Position = [l.Position(1)+0.09 l.Position(2)-0.01 l.Position(3) l.Position(4)];
%     end

pos = get(gca,'Position');

colormap(parula(length(phenmap)));
c=colorbar;
% % c.Limits = [1 1.9];
c.Ticks = [];
title(c,{'PC1'});
cpos = c.Position;
% cpos(1) = cpos(1)-0.1;
cpos(2) = cpos(2)+0.1;
cpos(4) = cpos(4)-0.2;
set(c,'Position',cpos);
set(gca,'Position',pos);

    
end
