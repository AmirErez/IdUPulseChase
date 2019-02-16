function data = PrepareForPairwise(expdata,DifferentiationPaths)
% Input: DifferentiationPaths has a cell for each tissue type
%        Each tissue had an ordered list of scatterslice coords for the linear path
% Goes over [xId,yId] in DifferentiationPaths, defining populations
% Averages across replicates and arranges time-series according to tabSamples
% prepares the 'data' struct which can be used for pairwise fitting

% Set timeseries parameters
data = struct;
data.time = unique(expdata.tabSamples.Time);
data.ntime = length(data.time);
Replicates = unique(expdata.tabSamples.Replicate);
Replicates = Replicates(~isnan(Replicates));

% Make populations
data.tabPopulations = table;
cnt = 1;
for tp=1:expdata.nType
    
    for ss=1:length(DifferentiationPaths{tp})
        xx = DifferentiationPaths{tp}(ss,1);
        yy = DifferentiationPaths{tp}(ss,2);
        
        popname = [char(expdata.Type(tp)) '_x' num2str(xx),'_y' num2str(yy)];
        tempdata.Population = categorical({popname});
        tempdata.Type = expdata.Type(tp);
        tempdata.Index = cnt;
        cnt = cnt+1;
        
        % Depricated
%         tempdata.Name = categorical({popname}); %redundant with Population
%         tempdata.Active = true;
%         tempdata.Min_children = 1;
%         tempdata.Max_children = 1;
%         tempdata.Hierarchy = 9;
%         tempdata.ColorR = 0.5;
%         tempdata.ColorG = 0.4;
%         tempdata.ColorB = 0.3;
%         tempdata.Proliferation_index = 0;
        
        tempdata.yPlotPos = yy/expdata.CoarseGrid.leny;
        tempdata.xPlotPos = xx/expdata.CoarseGrid.lenx;
        
        tempdata.xId = xx;
        tempdata.yId = yy;
        data.tabPopulations = [data.tabPopulations; struct2table(tempdata)];
    end
end

data.Populations = data.tabPopulations.Population;
data.nPopulations = length(data.Populations);
data.Type = data.tabPopulations.Type;
data.Types = unique(data.Type);
data.nTypes = length(data.Types);
data.CoarseGrid = expdata.CoarseGrid;

% Collect observables, filter time<0
Observables = {'CellCount','MFI_all','CellFrac'};
nObservables = length(Observables);
ObsStruct = struct;
for tp=1:expdata.nType
    for oo=1:nObservables
        tabname = [char(expdata.Type(tp)) '_' Observables{oo}];
        ObsStruct.(tabname) = table;
    end
end
for ss=1:height(expdata.tabSamples)
    multiWaitbar('PrepareForPairwise', ss/height(expdata.tabSamples));
    if(expdata.tabSamples.Time(ss)<0), continue; end
  
    tp = find(expdata.Type==expdata.tabSamples.Type(ss));
    tempstruct = struct;
    for oo=1:nObservables
        tabname = [char(expdata.Type(tp)) '_' Observables{oo}];
        tempstruct.Sample = expdata.tabSamples.Sample(ss);
        pops = data.tabPopulations.Population(data.tabPopulations.Type==expdata.Type(tp));
        
        ob = NaN(1,length(pops));
        for pp=1:length(pops)
            pop_ind = find(data.tabPopulations.Population==pops(pp));
%             warning('!!!Check binned samples ordered the same as populations');
            pop_eventind = find((expdata.Data{ss}.xId==data.tabPopulations.xId(pop_ind))...
                    .* (expdata.Data{ss}.yId==data.tabPopulations.yId(pop_ind)));
            if(isempty(pop_eventind))
                val = NaN;
            elseif(strcmp(Observables{oo},'CellCount'))
                val = nnz(pop_eventind);
            elseif(strcmp(Observables{oo},'CellFrac'))
                val = nnz(pop_eventind)/height(expdata.Data{ss});
            elseif(strcmp(Observables{oo},'PercentPositive'))
                val = 0;
            elseif(strcmp(Observables{oo},'gMFI_positive'))
                val = 0;
            elseif(strcmp(Observables{oo},'MFI_all'))
                val = nanmean(expdata.Data{ss}.(expdata.CoarseGrid.EdUChannel)(pop_eventind));
            else
                warning(['Unknown observable ' Observables{oo} ' in sample ' char(expdata.tabSamples.Sample(ss))]);
            end
            
            tempstruct.(char(pops(pp))) = val;            
        end
        ObsStruct.(tabname) = [ObsStruct.(tabname); struct2table(tempstruct)];
    end
    
    
end
 multiWaitbar('PrepareForPairwise','Close');
 
% Contruct time-series, average observables over replicates
for oo=1:nObservables
    
   data.(Observables{oo}) = NaN(data.ntime,length(Replicates),data.nPopulations);
   for tp=1:data.nTypes
       tabname = [char(expdata.Type(tp)) '_' Observables{oo}];
       tab = ObsStruct.(tabname);
       
       try
           tab = join(tab,expdata.tabSamples);
       catch ME
           error(['Cannot match ' tabname]);  
       end
       % get read populations
       pops = tab.Properties.VariableNames(2:(end-width(expdata.tabSamples)+1));
    
       ind_pops = zeros(length(pops),1);
       for pp=1:length(pops)
           try
               ind_pops(pp) = data.tabPopulations.Index(data.tabPopulations.Population == pops{pp});
           catch
               error(['Cannot find population ' pops{pp} '. Check Populations tab']);
           end
       end
          
       % Go sample by sample and read data into struct elements
       for ss=1:height(tab)
           tt = find(data.time==tab.Time(ss));
           rr = find(Replicates == tab.Replicate(ss));
           %     disp([rr tt]);
           if(isempty(tt) || isempty(rr) || isnan(rr) || rr==0)
               warning(['Skipping sample ' tab.Sample{ss}]);
           else
               data.(Observables{oo} )(tt,rr,ind_pops) = squeeze(tab{ss,2:(2+length(pops)-1)});
           end
       end
                         
   end
   
   % Get mean and err
   obs = Observables{oo};
   obs_err = [obs '_err'];
   data.(obs_err) = squeeze(nanstd(data.(obs),1,2))./sqrt(squeeze(nansum(~isnan(data.(obs)),2))-1);
   data.(obs) = squeeze(nanmean(data.(obs),2));
   if(strcmp(obs,'MFI_all')) % Remove background fluorescense 
       data.MFI_all = data.MFI_all - repmat(min(data.MFI_all),data.ntime,1); 
%        data.MFI_all = data.MFI_all - repmat(data.MFI_all(1,:),data.ntime,1);  
%        data.MFI_all = data.MFI_all - min(data.MFI_all(:));
%        data.MFI_all(data.MFI_all<0) = 1;
   end
   
end
 


%% Plot MFI
% phenmap = colormap(parula(length(expdata.CoarseGrid.y)));
% for tp=1:data.nTypes
%     figure;
%     hold on
%     set(gcf,'Name',char(data.Types(tp)));
%     xlabel('Time [hr]');
%     ylabel('MFI');
%     
%     usepopulations = data.tabPopulations.Index(data.tabPopulations.Type==data.Types(tp));
%     
%     for pp=1:length(usepopulations)
%        popind = usepopulations(pp);
%        errorbar(data.time,data.MFI_all(:,popind),data.MFI_all_err(:,popind),'.-','Color',phenmap(pp,:),'DisplayName',char(data.Populations(popind))); 
%     end
%     set(gca,'YScale','log')
% end