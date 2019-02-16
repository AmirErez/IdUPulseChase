function expdata = ImportTxtToExpdata(dirname)

%% Read data:
expdata = struct;


if(~isfield(expdata,'DataSuffix')), expdata.DataSuffix = '.csv'; end
if(~isfield(expdata,'Delimiter')), expdata.Delimiter = ','; end
% expdata.DataSuffix = '.txt';
% expdata.Delimiter = '\t';
expdata.dirname = dirname;
expdata.SamplesFile = [dirname filesep 'Samples.xlsx'];

expdata.tabSamples = readtable([expdata.SamplesFile]);
expdata.tabSamples.Type = categorical(expdata.tabSamples.Type);
if(any(strcmp('Tissue',fieldnames(expdata.tabSamples))))
    expdata.tabSamples.Tissue = categorical(expdata.tabSamples.Tissue);
end
expdata.Data = cell(height(expdata.tabSamples),1);
for ss=1:height(expdata.tabSamples)
    filename =[expdata.dirname filesep expdata.tabSamples.Sample{ss}] ;
    disp(['Reading ' filename]);
    multiWaitbar('Reading',ss/height(expdata.tabSamples));
    try
        tabData = readtable([filename expdata.DataSuffix],'Delimiter',expdata.Delimiter);
    catch
       disp(['Problem reading file ' filename expdata.DataSuffix]); 
    end
%     tabData = readtable([filename expdata.DataSuffix],'Delimiter',','); 
%     tabData = readtable([filename '.csv'],'Delimiter',','); 
    tabData{:,:} = tabData{:,:}+5*(1+2*(rand(size(tabData{:,:}))-0.5));
    expdata.Data{ss} = tabData;
    
    if(height(expdata.Data{ss})<100)
        warning(['Sample ' filename ' has less than 100 events !']);
    end
end
multiWaitbar('Reading','Close');

expdata.Data = expdata.Data(expdata.tabSamples.Time>=0); 
expdata.tabSamples = expdata.tabSamples(expdata.tabSamples.Time>=0,:);

