function expdata = ScatterSlice(expdata_,CoarseGrid)
% Slices expdata.Data according to CoarseGrid information, saves back in expdata.Data

expdata = expdata_;

for ss=1:height(expdata.tabSamples)          
    multiWaitbar('ScatterSlice',ss/height(expdata.tabSamples));
    val = log10(expdata.Data{ss}{:,{CoarseGrid.ychan,CoarseGrid.xchan}});
    yId = 1+sum( bsxfun( @ge, val(:,1), CoarseGrid.y(1:end-1)' ), 2 ) ;
    xId = 1+sum( bsxfun( @ge, val(:,2), CoarseGrid.x(1:end-1) ), 2 ) ;

    tempstruct = struct;
    tempstruct.xId = xId;
    tempstruct.yId = yId;
    expdata.Data{ss} = [expdata.Data{ss},struct2table(tempstruct)];
end

multiWaitbar('ScatterSlice','Close');
    

% Now plot
% Types = unique(expdata.tabSamples.Type);
% nTypes = length(Types);

expdata.ScatterSliceHists = cell(expdata.nType,1);

for tp=1:expdata.nType
    cellcounts  = zeros(expdata.CoarseGrid.leny,expdata.CoarseGrid.lenx);
   
    
    fn = find(expdata.tabSamples.Type == expdata.Type(tp));
    hsts = zeros(length(fn),length(CoarseGrid.y));
    for ss=1:length(fn)
        if(expdata.tabSamples.Type(ss)~=expdata.Type(tp)), continue; end
        tempcounts  = zeros(expdata.CoarseGrid.leny,expdata.CoarseGrid.lenx);
         
        for yy=1:expdata.CoarseGrid.leny
            for xx=1:expdata.CoarseGrid.lenx
                tempcounts(yy,xx)=nnz((expdata.Data{ss}.xId==xx).*(expdata.Data{ss}.yId==yy));
            end
        end
        cellcounts = cellcounts + tempcounts / sum(tempcounts(:));
    end
    
    expdata.CoarseGrid.ScatterSliceCellCounts{tp} = cellcounts;
    
   
    to_plot = false;
    if(to_plot)
        figure;
        imagesc(expdata.CoarseGrid.x, expdata.CoarseGrid.y,...
            log10(expdata.CoarseGrid.ScatterSliceCellCounts{tp}));
        set(gca,'YDir','Normal');
        colormap(parula(256));
        set(gca,'FontSize',18);
        title(char(expdata.Type(tp)));
        set(gcf,'Name',['ScatterSlice Density ' char(expdata.Type(tp))]);
        
        set(gca,'YTick',expdata.CoarseGrid.y(1:2:end));
%         set(gca,'YTickLabel',cellstr(['10^{' num2str(expdata.CoarseGrid.y(1:2:end)) '}']));
        set(gca,'XTick',expdata.CoarseGrid.x(1:2:end));
%         set(gca,'XTickLabel',cellstr(['10^{' num2str(expdata.CoarseGrid.x(1:2:end)) '}']));

        xlabel(CoarseGrid.xchanName);
        ylabel(CoarseGrid.ychanName);
%         ylim([0 max(meanhist)*1.3]);
    end
%     print(gcf,'-djpeg',[expdata.dirname filesep char(expdata.Type(tp)) 'Sliced'],'-r150');
end