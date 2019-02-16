function expdata = GateRelevantCells(expdata_, Gates)
expdata = expdata_;

fracs = zeros(height(expdata.tabSamples),1);

CumAbundances = cell(expdata.nType,1);


for tp=1:expdata.nType
    CumAbundances{tp} = zeros(256);
    X = zeros(256);
    Y = zeros(256);
    for ss=1:height(expdata.tabSamples)
        if(expdata.tabSamples.Type(ss)~=expdata.Type(tp)), continue; end
        % Get non-negative expression levels
        val = expdata.Data{ss}{:,{expdata.CoarseGrid.xchan,expdata.CoarseGrid.ychan}};
        NonNegative = find((val(:,1)>0).*(val(:,2)>0));
        val = log10(val);
        %     nonnegative = intersect(find(val(:,1)>0),find(val(:,2)>0));
        %     [y,~] = ind2sub(size(val),nonnegative);
        %     val = val(unique(y),:);
        
        %Plot
        to_plot = false;
        to_save = false;
        if(to_plot)
            ContourEvents(val(NonNegative,:),{expdata.CoarseGrid.xchanName,expdata.CoarseGrid.ychanName},char(expdata.tabSamples.Sample(ss)),to_save);
        end
        
        InGate = zeros(height(expdata.Data{ss}),length(Gates.Channels)); % 1 means in gate, 0 outside
        
        for cc=1:length(Gates.Channels)
            chan = Gates.Channels{cc};
            
            % Make sure we're inside the gate
            
            fn=find(  (expdata.Data{ss}.(chan)>=Gates.Ranges(cc,1)).*...
                (expdata.Data{ss}.(chan)<Gates.Ranges(cc,2))        );
            InGate(fn,cc) = 1;
            
        end
        
        % Find max density in gate
        EventsInGate = find(prod(InGate,2));
        
        ffn = intersect(NonNegative,EventsInGate); % Use on raw data
        
%         [~,density,X,Y]=kde2d(val(ffn,:));
        [~,density,X,Y]=kde2d(val(NonNegative,:));
        CumAbundances{tp} = CumAbundances{tp} + density;
        
        fracs(ss) = length(ffn)/height(expdata.Data{ss});
        expdata.Data{ss} = expdata.Data{ss}(ffn,:);
        disp(char(expdata.tabSamples.Sample(ss)));
        disp(['   Using %' num2str(fracs(ss)*100) ' of input.']);
    end
    
    CumAbundances{tp} =  CumAbundances{tp}/sum( CumAbundances{tp}(:));    
    density = CumAbundances{tp};
    density(density<0) = 0;
    density = density / max(density(:));
    density(density<10^-4) = NaN;
    [mx mxind]= max(density,[],2);
    subplot(1,expdata.nType,tp);
    set(gca,'FontSize',10);
    imagesc(X(1,:),Y(:,1),log10(density));
    set(gca,'YDir','Normal');
    c=colorbar('Location','West');
    xlabel(expdata.CoarseGrid.xchanName);
    ylabel(expdata.CoarseGrid.ychanName);
%     s=surf(X(1,:),Y(:,1),log10(density));
%     hold on
%     s.EdgeAlpha = 0.15;
%     title(char(expdata.Type(tp)));  
%     plot3(X(1,mxind),Y(:,1),density(:,mxind),'ok','MarkerSize',6,'MarkerFaceColor','w');
%     plot(X(1,mxind),Y(:,1),'xk','MarkerSize',6,'MarkerFaceColor','w');
        
end
  disp(['Summary: Using %' num2str(nanmean(fracs)*100,3) char(177) num2str(nanstd(fracs*100),3) ' of input.']);
