function fitout = fitPairParseArgs(time, source_positive_fraction,dest_positive_fraction, dest_err, varargin)
% Fits dynamics of source to destination for nSteps steps.
% Finds the d's and kappa's for each of the steps

p = inputParser;
addRequired(p,'time',@isnumeric);
addRequired(p,'source_positive_fraction',@isnumeric);
addRequired(p,'dest_positive_fraction',@isnumeric);
addRequired(p,'dest_err',@isnumeric);
addOptional(p,'nSteps',1,@isnumeric);
addOptional(p,'WithDelay',false,@islogical);
addOptional(p,'BruteForce',false,@islogical); % for nSteps=1, optimizes p,d by brute force
addOptional(p,'CostLandscapeSavefile',false,@isstr); % if BruteForce, saves to filename

parse(p,time, source_positive_fraction,dest_positive_fraction, dest_err,varargin{:})

time = p.Results.time;
source_positive_fraction = p.Results.source_positive_fraction;
dest_positive_fraction = p.Results.dest_positive_fraction;
dest_err = p.Results.dest_err;
nSteps = p.Results.nSteps;
WithDelay = p.Results.WithDelay;

if(p.Results.BruteForce && p.Results.nSteps>1)
    throw(MException('fitPairParseArgs:incompatible','Cannot use BruteForce with nSteps>1'));
end
if(p.Results.BruteForce && p.Results.WithDelay==true)
    throw(MException('fitPairParseArgs:incompatible','Cannot use BruteForce with WithDelay'));
end
fitout = struct;
fitin = struct;



% Record inputs
fitin.with_delay = WithDelay;
fitin.time = time;
fitin.source = source_positive_fraction;
fitin.dest = dest_positive_fraction;
fitin.dest_err = dest_err;
fitin.nSteps = nSteps;
fitin.BruteForce = p.Results.BruteForce;

% Record fit out
if(~fitin.BruteForce)
    fitout.fitopts = optimset('MaxIter',500000,'MaxFunEvals',50000,'TolFun',10^-10);
end
fitout.fitin = fitin;
fitout.kappas = NaN(nSteps,1);
fitout.ds = NaN(nSteps,nSteps);
fitout.delays = NaN(nSteps,nSteps);
fitout.SSE = NaN(nSteps,1);
fitout.SSEnormbyErr = NaN(nSteps,1);
fitout.exitflags = NaN(nSteps,1);
fitout.Lags = 4:5;
fittout.LBQ_pValue = NaN(nSteps,length(fitout.Lags));
fitout.residuals_normby_err = NaN(length(time),nSteps);
fitout.residuals_normby_err = NaN(length(time),nSteps);
fitout.unnorm_residuals = NaN(length(time),nSteps);
fitout.Jacobian = NaN; % Jacobian
fitout.confint = NaN; % 95% confidence intervals

% Remove NaNs from analysis. some problems not fully debugged
if(~isempty(find(isnan(source_positive_fraction),1)) || ~isempty(find(isnan(dest_positive_fraction),1)))
    warning('fitPairWithWithoutDelay: Found NaN in input. Aborting.');
    return;
end
% notnans = boolean((~isnan(source_positive_fraction)).*(~isnan(dest_positive_fraction)));
% time = time(notnans);
% source_positive_fraction = source_positive_fraction(notnans);
% dest_positive_fraction = dest_positive_fraction(notnans);
% dest_err = dest_err(notnans);

% Cost function as sum of squared errors (fitout.SSE):
if(fitin.with_delay==false)
    model = @(args,steps) ((calcSourceDestN(time,source_positive_fraction,dest_positive_fraction(1),args(1),args(2:end),steps)-dest_positive_fraction)...
    ./dest_err);
else
   model = @(args,steps) ((calcSourceDestNWithDelay(time,source_positive_fraction,dest_positive_fraction(1),args(1),args(2:end-1),args(end),steps)-dest_positive_fraction)...
    ./dest_err); 
end

modelvalonestep = @(args,src) calcSourceDestN(time,src,dest_positive_fraction(1),args(1),args(2:end),length(args(2:end)));

for ss=1:nSteps
    if(~fitin.with_delay)
        if(~p.Results.BruteForce)
            args = zeros(ss+1,1);
            lbargs = zeros(length(args),1);
%             ubargs = 10*ones(length(args),1);
            ubargs = 1*ones(length(args),1);
            [estimates,fitout.SSE(ss),residuals,fitout.exitflags(ss),out,lambda,Jacobian] = lsqnonlin(@(args)model(args,ss),args,...
                lbargs,ubargs,fitout.fitopts);
            
            
            fitout.Jacobian = Jacobian;
            fitout.confint = nlparci(estimates,residuals,'Jacobian',Jacobian,'alpha',0.5);
            fitout.kappas(ss) = estimates(1);
            fitout.ds(1:ss,ss) = estimates(2:(ss+1));
            fitout.kappas_err(ss) = fitout.confint(1,2)-estimates(1);
            fitout.ds_err(1:ss,ss) = fitout.confint(2:end,2)-estimates(2:(ss+1));
            
            [ypred,delta] = nlpredci(modelvalonestep,source_positive_fraction,estimates,residuals,'Jacobian',full(Jacobian));
            fitout.ypred = ypred;
            fitout.delta = delta;
            DOF = fitout.Lags - 2;
        else
            % Brute force method ! Which means check all p,d values in [0 1] range
            if(nSteps>1)
               throw(MException('fitPairParseArgs:incompatible','Cannot use BruteForce with nSteps>1'));
            end
            if(WithDelay)
                throw(MException('fitPairParseArgs:incompatible','Cannot use BruteForce with delay'));
            end
            allks = [0:0.01:0.5];
            allds = [0:0.01:0.5];
            
            normSSE = zeros(length(allks),length(allds));
            vals = cell(size(normSSE));
            for kk=1:length(allks)
                multiWaitbar('Optimizing',kk/(length(allks)));
                for dd=1:length(allds)
                    vals{kk,dd} = calcSourceDestN(time,source_positive_fraction,dest_positive_fraction(1),allks(kk),allds(dd),1);
                    
                    residuals = dest_positive_fraction - vals{kk,dd};
                    normSSE(kk,dd) = sum((residuals.^2) ./ dest_err.^2);                    
                end
            end
            [mn,mnind] = min(normSSE(:));
            residuals = (dest_positive_fraction - vals{mnind})./dest_err; % "residuals" always norm by error
            [yind,xind] = ind2sub(size(normSSE),mnind);
            fitout.kappas = allks(yind);
            fitout.kappas_err = 0;
            fitout.ds = allds(xind);
            fitout.ds_err = 0;
            fitout.delays = NaN;
            fitout.ypred = vals{mnind};
            fitout.delta = 0;

            
            DOF = fitout.Lags - 2;
%             figure; 
%             imagesc(allks,allds,log10(normSSE));
%             xlabel('k');
%             ylabel('d');            
            multiWaitbar('Optimizing','Close');
            
        end
    else
        args = zeros(ss+2,1);
        ubargs =10*ones(length(args),1);
        ubargs(3) = max(time)-1;
        lbargs = zeros(length(args),1);
        
        [estimates,fitout.SSE(ss),residuals,fitout.exitflags(ss),out,lambda,Jacobian] = lsqnonlin(@(args)model(args,ss),args,...
            lbargs,ubargs,fitout.fitopts);
        fitout.Jacobian = Jacobian;
        fitout.confint = nlparci(estimates,residuals,'Jacobian',J);
        fitout.kappas(ss) = estimates(1);
        fitout.ds(1:ss,ss) = estimates(2:(ss+1));
        fitout.delays(ss) = estimates(end);
        
        DOF = fitout.Lags - 3;
    end
    fitout.residuals_normby_err(:,ss) = residuals;
    fitout.unnorm_residuals(:,ss) = residuals.*dest_err;

    fitout.SSEnormbyErr(ss) = sum((fitout.residuals_normby_err(:,ss).^2));
    fitout.SSEunnorm(ss) = sum((fitout.unnorm_residuals(:,ss).^2));
    
    % Resample residuals with equi-spaced time-points for LBQ-test. Otherwise 'lag' is poorly defined
    resampled_residuals = interp1(time,fitout.residuals_normby_err(:,ss),linspace(time(1),time(end),length(time)));
    [~,pValue,~,~] = lbqtest(resampled_residuals,'lags',fitout.Lags,'DOF',DOF,'Alpha',0.001);
    fitout.LBQ_pValue(ss,:) = pValue;
end



if(~isempty(find(fitout.exitflags<0,1)))
    error('Bad optimization !');
end

nTime = length(time);
% K = [1:nSteps]' + 2; % Number of params: 1 for k, 1 for each step (d), and 1 for noise variance
K = 3*(1:nSteps)' ; % Number of params: 1 for k, 1 for each step (d), and 1 for noise variance

% minusLogLikelihood = 1/2*nTime*log(fitout.SSE/nTime); % I don't get this one
minusLogLikelihood = fitout.SSEnormbyErr(ss)/2;
% To calculate corrected AIC:
AIC = 2*minusLogLikelihood+2*K;
fitout.AICc =  AIC + 2*K.*(K+1)./(nTime-K-1);

% fitout.AICc = minusLogLikelihood;
