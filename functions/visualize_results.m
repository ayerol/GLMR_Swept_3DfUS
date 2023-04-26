function visualize_results(U,V,t_axis,stim_times,Nz,Nx,Ny)

cmap = lines(7);
rank = size(U,2);

for r = 1:rank

    if ~isempty(V)
        figure;
        if r == 1
            plot(t_axis,V(:,r),'LineWidth',5,'Color',cmap(r,:));
        else
            p = plot(t_axis,V(:,r),'o','Color',cmap(r,:));
            p.MarkerFaceColor = cmap(2,:);
            p.MarkerSize = 7;
        end
        set(gca,'FontSize',26); xlabel('Time (s)'); xticks(0:30:180);
        xlim([0 180]);
        set(gca,'fontname','times')
        % ylabel('Amplitude')
    
        hold on;
        if r == 1
            ymax = 1.03*max(V(:,1)); ymin = 0;
        else
            ymax = 4.9; ymin = -4;
        end
    
        for i = 1:6
            p3 = fill([stim_times(i*2-1) stim_times(i*2-1) ...
                stim_times(i*2) stim_times(i*2)],...
                [ymin ymax ymax ymin],'r','HandleVisibility','off');
            set(p3,'facealpha',.1);
            set(p3,'EdgeColor','none');
        end
        ylim([ymin ymax]);
    end

    figure;
    h = slice(permute(reshape(U(:,r),Nz,Nx,Ny),[1 3 2]),1:Ny,[],[]);
    set(h,'EdgeColor','none',...
        'FaceColor','interp',...
        'FaceAlpha','interp')
    alpha('color')
    xlabel('Slice index','FontSize',25)
    zlabel('Depth','FontSize',25)
    ylabel('Width','FontSize',25)
    set(gca,'FontSize',22);
    alphamap('increase',.1) %.05 in paper figures
    colormap parula;
    caxis([-.005 max(U(:,r))]);

    if r == 1
        cb = colorbar;
        cb.Location = 'north';
        cb.Position = [0.1443 0.84 0.76 0.0381];
        cbtitle = get(cb,'Title');
        set(cbtitle,'String','Estimated Task-Response Strength');
        set(cbtitle,'Position',[212.8 68 0]);
    else
        cb = colorbar;
        cb.Location = 'north';
        cb.Position = [0.1443 0.79 0.76 0.0381];
        cbtitle = get(cb,'Title');
        set(cbtitle,'String','Estimated Backgr. Activity Strength');
        set(cbtitle,'Position',[212.8 64 0]);
    end

    ax = gca;
    ax.Position = ax.Position - [0 0 0 .15];
    zticks([0 10 20]); xticks([0 10 20]); yticks([0 10 20]);
    set(cbtitle,'FontSize',24);
    set(gca,'fontname','times')

end

end