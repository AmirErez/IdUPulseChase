mu = [2,3];
sigma = [1,1.5;1.5,5];
% rng default  % For reproducibility
r = mvnrnd(mu,sigma,1000);

figure
hold on
plot(r(:,1),r(:,2),'+');

c = cov(r);
[vec,val] = eig(c);

[~,mx] = max(diag(val));
largestvec = vec(:,mx);
angl = atan(largestvec(2)/largestvec(1));

mn = mean(r);
len = 4;
plot(mn(1)+[0,cos(angl)*len],mn(2)+[0,sin(angl)*len],'-k','LineWidth',2);
plot(mn(1)-[0,cos(angl)*len],mn(2)-[0,sin(angl)*len],'-k','LineWidth',2);

