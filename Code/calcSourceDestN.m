function dest = calcSourceDestN(time, source_positive_fraction,dest_positive_initial, kappa,d,N)
% Calculate dest for source N steps removed. N=1,2,3 for now
% According to recursion formula outlined in the code
% d is size Nx1 is disappearence rate of each generation (0..N).
% kappa is size Nx1 which is composed of the "differentiation rates" 
%       so that [k_10 k_20 k_30] are the one/two/three step rates from source to destination
% dest is size time x N
% dest_positive_initial is the initial value of the destination population (at first timepoint)
%                       this allows to better fit starting from time >0

nTime = length(time);
dest = zeros(nTime,1);
d = d(1:N);

if(N<1 || N>5), error('Bad N'); end


% Get Delta which is the convolution integrals with all the ds
Delta = zeros(nTime,N);
repdt = repmat(d',nTime,1).*repmat(time,1,N);
Delta = exp(-repdt).*cumtrapz(time,repmat(source_positive_fraction,1,N).*exp(repdt));
    
if(N==1) %always calculate
    %     F_1(t) = \kappa_10 \int_0^t dtprime F_0(t_0)*exp(d*tPrimeMinust)
    dest = kappa*Delta(:,1);
%     dest_direct = kappa*exp(-d(1)*time).*cumtrapz(time,source_positive_fraction.*exp(d(1)*time));
    
elseif(N==2)
    %     F_2(t) =
% Calculate by a single integral "Delta" because of inverting the double-integral order 
% and integrating one integral:
    %     dest = kappa/(d(2)-d(1))*((Delta(:,1)-Delta(:,2)));
    
% Calculate directly by integrating the 1-step to give the 2-step result
    dest1 = Delta(:,1);
    dest2 = exp(-d(2)*time).*cumtrapz(time,dest1.*exp(d(2)*time));
    dest = kappa*dest2;
    
% Calculate by converting (Delta1 - Delta2) / (d2-d1) to sinhc
% For some reason this is numerically unstable 
%     su = (d(2)+d(1))/2;
%     di = (d(2)-d(1))/2;
%     for tt=2:length(time)
%         tminustprime = time(tt)-time(1:tt);
%         dest(tt) = trapz(time(1:tt),exp(-su*tminustprime).*sinhc(di*tminustprime).*tminustprime);
%     end
%     dest = dest*kappa;

elseif(N==3)
    %     F_3(t) =
    % Calculate from the Delta integral:
%     dest = kappa*abs(...
%         (d(3)-d(2))*Delta(:,1) + (d(1)-d(3))*Delta(:,2) + (d(2)-d(1))*Delta(:,1) );
   
    % Calculate directly by integrating twice the 1-step solution
    dest2 = exp(-d(2)*time).*cumtrapz(time,Delta(:,1).*exp(d(2)*time));
    dest3 = exp(-d(3)*time).*cumtrapz(time,dest2.*exp(d(3)*time));
    dest = kappa*dest3;
    
elseif(N==4)
    % Using Delta i:
%     dest = kappa/(d(2)-d(1))*...
%         abs( ( (Delta(:,1)-Delta(:,4))/(d(4)-d(1)) - (Delta(:,3)-Delta(:,4))/(d(4)-d(3))  /(d(3)-d(1))) ...
%          -( (Delta(:,2)-Delta(:,4))/(d(4)-d(2)) - (Delta(:,3)-Delta(:,4))/(d(4)-d(3))  /(d(3)-d(2)))  );
    dest2 = exp(-d(2)*time).*cumtrapz(time,Delta(:,1).*exp(d(2)*time));
    dest3 = exp(-d(3)*time).*cumtrapz(time,dest2.*exp(d(3)*time));
    dest4 = exp(-d(4)*time).*cumtrapz(time,dest3.*exp(d(4)*time));
    dest = kappa*dest4;
elseif(N==5)
    dest2 = exp(-d(2)*time).*cumtrapz(time,Delta(:,1).*exp(d(2)*time));
    dest3 = exp(-d(3)*time).*cumtrapz(time,dest2.*exp(d(3)*time));
    dest4 = exp(-d(4)*time).*cumtrapz(time,dest3.*exp(d(4)*time));
    dest5 = exp(-d(5)*time).*cumtrapz(time,dest4.*exp(d(5)*time));
    dest = kappa*dest5;
    
end
    
  dest = dest + dest_positive_initial*exp(-(time-time(1))*d(end));    

end
    
