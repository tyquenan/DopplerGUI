function pca_ripple(nrem,eigen,curpath,dx,dy,mx,scores,EXP)
    EV_nrem = find([nrem.Value] == 1);
    load(string(curpath)+'Time_Reference.mat');
    load(string(curpath)+'Time_Groups.mat');
    save(string(curpath)+'EVnrem.mat', 'EV_nrem');
    t = time_ref.Y;
    ripp = table2array(readtable(string(curpath) + "Events"+ filesep + "Ripples-Abs-All.csv"));
    t_events = ripp(:,1) - t(1);
    t=t-t(1);
    n_events = length(t_events);
    % Computing event averages and fUS averages
    t_before = -1;          % time window before (seconds)
    t_after = 5;            % time window after (seconds)
    sampling_fus = 10;      % Hz
    %     % comment when done
    %     sampling_fus = 5;      % Hz
    t_bins_fus  = (t_before:1/sampling_fus:t_after)';
    n_time = length(t_bins_fus);
    % Interpolate fUS
    Xq_evt_fus = [];
    for i =1:length(t_events)
        Xq_evt_fus = [Xq_evt_fus;t_events(i)+t_bins_fus];
    end
    f = figure;
    f0=uitabgroup(f,'Units','Normalized','Position',[0 0 1 1]);
    dvrem=length(EV_nrem);
    f2=uitab(f0,'Units','Normalized','Title',"film");
    f=uitab(f0,'Units','Normalized','Title',"component");        
    XX = zeros(dvrem,n_time,dx*dy);
    for i = 1:dvrem
        X=zeros(n_events,n_time,dx*dy);
        ax = axes(f,'Position',[0.005+mod(i-1,4)/4,0.75-0.25*floor((i-1)/4),1/(4.2),0.24]);
        Y2q_evt = (interp1(t,scores(EV_nrem(i),:),Xq_evt_fus));
        mx_evt = (interp1(t,mx,Xq_evt_fus));
        mx_evt = reshape(mx_evt,[length(t_bins_fus) length(t_events)]);
        Y2q_evt = reshape(Y2q_evt,[length(t_bins_fus) length(t_events)]);
        imagesc('XData',t_bins_fus,'YData',1:length(t_events),'CData',Y2q_evt',"Parent",ax);
        axis off tight;
        for rip = 1:n_events
            X(rip,:,:)= Y2q_evt(:,rip)*eigen(EV_nrem(i),:)+mx_evt(:,rip);
        end
        XX(i,:,:) = nanmean(X,1);
    end
    video = VideoWriter(curpath+"ripples.avi"); %Create a video object
    video.FrameRate = 10;
    open(video); % Open video source - restricts the use of video for your program
    CLs=zeros(dvrem,2);
    for m=1:n_time
        axtxt= axes(f2,'Position',[.2,.95,.7,.05]);
        axis off
        txt = ['Ripples' EXP ', time from event peak:' num2str((m-1-(n_time/6))*6/n_time) 's'];
        text(axtxt,0,0,txt);
        for i=1:dvrem
            ax=axes(f2,'Position',[0.005+mod(i-1,4)/4,0.65-0.25*floor((i-1)/4),1/(4.2),0.24]);
            ax.XLim =[0 dy];
            ax.YLim = [0 dx];
            regg=col2im(squeeze(XX(i,m,:)),[1 1],[dx dy]);
            imagesc(ax,regg);
            ax.XTickLabel =''; 
            ax.YTickLabel =''; 
            axis off

            if m==1
                CLs(i,:) =ax.CLim;
            else
                ax.CLim = CLs(i,:);
            end
            colorbar;
        end
        drawnow;
        vidFrame = getframe(gcf);
        writeVideo(video,vidFrame); 
        cla(axtxt);
    end
save(curpath +"ripple.mat","XX","Y2q_evt");
end