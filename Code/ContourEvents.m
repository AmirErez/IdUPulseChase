function CountourEvents(val,Channels,Title,to_save)

figure;
set(gca,'FontSize',18);
hold on
grid('on');
[~,density,X,Y]=kde2d(val);
%         imagesc(X(1,:),Y(:,1),log(density+1));
contour(X(1,:),Y(:,1),log(density+1),100);
%             s=surf(X(1,:),Y(:,1),log(density+1));
%             s.EdgeColor = 'w';
%             s.EdgeAlpha = 0.15;
colormap(parula(256));
colorbar;
title(Title,'Interpreter','None');

set(gca,'YDir','Normal');
xlabel(Channels{1},'Interpreter','None');
ylabel(Channels{2},'Interpreter','None');
if(to_save)
    filename = ['Density_' Title ];
    print(gcf,'-djpeg',[expdata.dirname filesep filename '.jpg'],'-r150');
    close(gcf);
end