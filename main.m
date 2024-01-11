%% Begin
clearvars
close all force
clc
dv = 20;
alt=1;
rec=1;
files = uipickfiles('REFilter','.mat'); % cell {}
nb_files = length(files);
for ii = 1:nb_files
    file_path{ii} = files{ii}(1:max(strfind(files{ii},filesep)));
end

A= split(files,filesep);
filename = squeeze({A{end-nb_files+1:end}});
if length(files)==1
    EXP = A{end-1};
else
    EXP = {A{:,:,end-1}};
end

for cur=1:length(files)
    addpath(file_path{cur});
    if filename{cur} == "analyse.mat"
        load(string(file_path{cur})+string(filename{cur}));
        load(string(file_path{cur})+"Doppler.mat");
        film = Doppler_film;

    elseif filename{cur} == "EV_done.mat"
        load(string(file_path{cur})+string(filename{cur}));
        load(string(file_path{cur})+"Doppler.mat");
        film = Doppler_film;
        pca_reconstr(film, eigen,scores,rem,nrem,wake,noise,mx,phases,file_path{cur});
    else
        load(string(file_path{cur})+string(filename{cur}));
        film = Doppler_film;
        [eigen,scores,mx]= pca_t(film,dv);
        save(string(file_path{cur})+'analyse.mat','eigen','scores','mx');
    end
    EV_liste=[];
    if alt
        [eigen,scores,mx] = pca_t_cov_alt(film,eigen,scores,mx,dv);
    else
        [eigen,scores,mx] = pca_t_cov(film,eigen,scores,mx,dv);
    end
    dt = length(film(1,1,:));
    dx = length(film(:,1,1));
    dy = length(film(1,:,1));
    load('Time_Reference.mat');
    t = time_ref.Y;
    t=t-t(1);

    for i=1:dv
        figures(i)=uifigure('Units','normalized','Name',"EV "+num2str(i),"IntegerHandle","on");
        ax = uiaxes(figures(i),'Units','normalized','Position',[0.01,0.95-0.4*16/9,0.4*dy/dx,0.4*16/9]);
        ax.XLim =[0 dy];
        ax.YLim = [0 dx];
        regg=col2im(squeeze(eigen(i,:)),[1 1],[dx dy]);
        imagesc(ax,regg);
        %colorbar;
        ax=uiaxes(figures(i),'Units','normalized','Position',[0.03+0.4*dy/dx,0.7,0.44,0.2]);
        plot(ax,t,abs(scores(i,:)));
        avg = nanmean(scores(i,:));
        sig = nanstd(scores(i,:));
        miny =avg-2*sig;
        maxy = avg+2*sig;
        ax.YLim=[miny maxy];
        ax.XLim(2)=t(end);
        contrib = zeros(6,1);
        labels = cell(6,1);
        load('Time_Groups.mat');
        colors = ["blue","cyan","yellow","red","red","red"];
        phases =cell(dt,1);
        for ii=1:4
            siz= size(TimeGroups_S(ii).TimeTags_images);
            labels{ii}=string(TimeGroups_name(ii));
            for j=1:siz(1)
                phases((TimeGroups_S(ii).TimeTags_images(j,1)):(TimeGroups_S(ii).TimeTags_images(j,2)))=TimeGroups_name(ii);
                x=[(t(TimeGroups_S(ii).TimeTags_images(j,1))) (t(TimeGroups_S(ii).TimeTags_images(j,2))) (t(TimeGroups_S(ii).TimeTags_images(j,2))) (t(TimeGroups_S(ii).TimeTags_images(j,1)))];
                y=[miny miny maxy maxy];
                patch(ax,x,y,colors(ii),'FaceAlpha',0.5,'EdgeColor','none');
                contrib(ii) = contrib(ii)+nansum(abs(scores(i,TimeGroups_S(ii).TimeTags_images(j,1):TimeGroups_S(ii).TimeTags_images(j,2))));
            end
            dur = TimeGroups_duration{ii,1};
            contrib(ii)= contrib(ii)/str2double(TimeGroups_frames(ii));

        end
        ax=uiaxes(figures(i),'Units','normalized','Position',[0.05+0.4*dy/dx,0,0.4,0.6]);
        p=pie(ax,contrib);
        for iii=1:4
            p(2*iii - 1).FaceColor=colors(iii);
        end
        btn_r(i) = uibutton(figures(i),"State","Text","REM","Position",[10 80 80 20]);
        btn_nr(i) = uibutton(figures(i),"State","Text","NREM","Position",[10 55 80 20]);
        btn_w(i) = uibutton(figures(i),"State","Text","Wake","Position",[10 30 80 20]);
        btn_noise(i) = uibutton(figures(i),"State","Text","Noise","Position",[10 5 80 20]);

        btn_val(i) =uibutton(figures(i),"Text","Valider","Position",[100 20 80 20],"ButtonPushedFcn", @(btn,event) valider(figures(i)));
    end
    fig = uifigure;
    uibutton(fig,"Text","Reconstruction","Position",[10 10 100 100],"ButtonPushedFcn",@(btn,event) pca_reconstr(film, eigen,scores,btn_r,btn_nr,btn_w,btn_noise,mx,phases,file_path{cur}));
end
