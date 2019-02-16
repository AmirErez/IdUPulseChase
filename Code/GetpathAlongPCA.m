function [temppath,AlongPCA] = GetpathAlongPCA(expdata,s,densityout)
   
   % Get main differentiation path:
   linex = densityout.mean(1)+cos(densityout.angl)*s;
   liney = densityout.mean(2)+sin(densityout.angl)*s;   
   [coarsex,coarsey,temppath] = GetpathAlongPCAHelper(expdata,linex,liney);
   
   % Get path points normalized along PCA
   AlongPCA = struct;
   AlongPCA.sPCA = interp1(linex,s,expdata.CoarseGrid.x);
   AlongPCA.sPCA  = AlongPCA.sPCA(~isnan(AlongPCA.sPCA));
   AlongPCA.sPCA  = AlongPCA.sPCA-min(AlongPCA.sPCA);
   AlongPCA.sPCA  = AlongPCA.sPCA/max(AlongPCA.sPCA);
   AlongPCA.coarsex = coarsex;
   AlongPCA.coarsey = coarsey;
   
   % Move one square perp to PC1 (ie. along PC2)   
%    shiftedx = linex + cos(densityout.angl-pi/2)*mean(diff(AlongPCA.sPCA));
%    shiftedy = liney + sin(densityout.angl-pi/2)*mean(diff(AlongPCA.sPCA));
%    [~,~,AlongPCA.shift1] = GetpathAlongPCAHelper(expdata,shiftedx,shiftedy);
%    shiftedx = linex - cos(densityout.angl-pi/2)*mean(diff(AlongPCA.sPCA));
%    shiftedy = liney - sin(densityout.angl-pi/2)*mean(diff(AlongPCA.sPCA));
%    [~,~,AlongPCA.shift2] = GetpathAlongPCAHelper(expdata,shiftedx,shiftedy);
    
%    grid2ind = reshape(1:(expdata.CoarseGrid.lenx*expdata.CoarseGrid.leny),...
%        expdata.CoarseGrid.leny,expdata.CoarseGrid.lenx);
%    
%    valspath = sub2ind([expdata.CoarseGrid.leny,expdata.CoarseGrid.lenx],temppath(:,1),temppath(:,2));
%    valshift1 = sub2ind([expdata.CoarseGrid.leny,expdata.CoarseGrid.lenx],shift1(:,1),shift1(:,2));  
   
   % Helper function
   function [coarsex,coarsey,temppath] = GetpathAlongPCAHelper(expdata,linex,liney)
       
%        meshx = repmat(linex,length(liney),1);
%        meshy = repmat(liney',1,length(linex));
%        
%        qcx = repmat(expdata.CoarseGrid.x,length(expdata.CoarseGrid.y),1);
%        qcy = repmat(expdata.CoarseGrid.y,1,length(expdata.CoarseGrid.x));
%     
%        tcx = repmat(expdata.CoarseGrid.x',1,length(expdata.CoarseGrid.y));    
%        tcy = repmat(expdata.CoarseGrid.y',length(expdata.CoarseGrid.x),1);
%      
%        coarsex = interp2(meshx,meshy,meshx,qcx,qcy);       
%        
%        coarsey = interp2(meshx,meshy,meshy,tcy,tcx)';
%        
%        
%        for jj=1:size(coarsex,2)
%            f = find((~isnan(coarsex(:,jj)).*~isnan(coarsey(:,jj))));
%            if(~isempty(f)),break; end
%        end
%        if(isempty(f)), warning('Error!'); end
%        coarsex = coarsex(:,jj);
%        coarsey = coarsey(:,jj);
       
       coarsex1 = interp1(linex,linex,expdata.CoarseGrid.x);
       coarsex2 = interp1(liney,linex,expdata.CoarseGrid.y);
       coarsex = coarsex1;
       if(length(~isnan(coarsex2))>length(~isnan(coarsex1)))
          coarsex = coarsex2;
       end
       coarsex = coarsex(~isnan(coarsex));
       
       coarsey1 = interp1(linex,liney,expdata.CoarseGrid.x);
       coarsey2 = interp1(liney,liney,expdata.CoarseGrid.y);
       coarsey = coarsey1;
       if(length(~isnan(coarsey2))>length(~isnan(coarsey1)))
          coarsey = coarsey2;
       end
       coarsey = coarsey(~isnan(coarsey));
       
       temppath = zeros(length(coarsex),2);
       for cc=1:length(coarsex)
           [~,xind] = min(abs(expdata.CoarseGrid.x-coarsex(cc)));
           [~,yind] = min(abs(expdata.CoarseGrid.y-coarsey(cc)));
           temppath(cc,:) = [xind,yind];
       end
   end

end