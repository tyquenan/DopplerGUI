%% Begin
clearvars
close all force
clc
dv = 20; %Number of PC
rec=1;
files = uipickfiles('REFilter','.mat'); % open the GUI to select files

%Supported formats: Doppler.mat, analyse.mat (PCA done), EV_done.mat (EV already
%selected)
%% Organize files
nb_files = length(files);
for ii = 1:nb_files
    file_path{ii} = files{ii}(1:max(strfind(files{ii},filesep)));
end
exp = zeros(nb_files,dv);
A= split(files,filesep);
filename = squeeze({A{end-nb_files+1:end}});
if length(files)==1 %check if only one file is selected
    EXP = A{end-1};
else
    EXP = {A{:,:,end-1}};
end
%% Mode selection (ALL, Only REM or Only NREM)
ftemp = uifigure;
ball = uibutton(ftemp,"State","Text","ALL","Position",[10 100 100 100]);
brem = uibutton(ftemp,"State","Text","REM","Position",[120 100 100 100]);
bnrem = uibutton(ftemp,"State","Text","NREM","Position",[230 100 100 100]);
while ball.Value == 0 && brem.Value == 0 && bnrem.Value == 0
    pause(0.1);
end
ball = ball.Value;
brem = brem.Value;
bnrem = bnrem.Value;
close ftemp
%% Main Loop
for cur=1:length(files)
    addpath(file_path{cur});
    load(string(file_path{cur})+'Time_Groups.mat');
    if filename{cur} == "analyse.mat" %PCA Already done
        load(string(file_path{cur})+string(filename{cur}));
        film = load(string(file_path{cur})+"Doppler.mat").Doppler_film;
        %% Apply selected mode
        if bnrem 
            for ii=1:4         
                if ii ~= 3
                    siz= size(TimeGroups_S(ii).TimeTags_images);
                    for j = 1:siz(1)
                        film(:,:,TimeGroups_S(ii).TimeTags_images(j,1):TimeGroups_S(ii).TimeTags_images(j,2)) = NaN;
                    end
                end
            end
        end
    elseif filename{cur} == "EV_done.mat" %PC already selected (careful on which mode !)
        load(string(file_path{cur})+string(filename{cur}));
        film = load(string(file_path{cur})+"Doppler.mat").Doppler_film;
        if bnrem
            for ii=1:4
                if ii ~= 3
                    siz= size(TimeGroups_S(ii).TimeTags_images);
                    for j = 1:siz(1)
                        film(:,:,TimeGroups_S(ii).TimeTags_images(j,1):TimeGroups_S(ii).TimeTags_images(j,2)) = NaN;
                    end
                end
            end
        end
        %Here only the reconstruction is displayed
        pca_reconstr2(eigen,scores,EV_rem,EV_nrem,EV_wake,EV_noiseless,mx,phases,file_path{cur});
    else
        %%Doppler.mat -> Normal mode
        load(string(file_path{cur})+string(filename{cur}));
        film = load(string(file_path{cur})+"Doppler.mat").Doppler_film;
        if bnrem
            for ii=1:4         
                if ii ~= 3
                    siz= size(TimeGroups_S(ii).TimeTags_images);
                    for j = 1:siz(1)
                        film(:,:,TimeGroups_S(ii).TimeTags_images(j,1):TimeGroups_S(ii).TimeTags_images(j,2)) = NaN;
                    end
                end
            end
        end
        [eigen,scores,latent,explained,mx]= pca_t(film,dv);%apply PCA
        save(string(file_path{cur})+'analyse.mat','eigen','scores','latent','explained','mx');
    end
    %% The user selects the PCs
    EV_liste=[];
    [eigen,scores,mx] = pca_t_cov_alt(film,eigen,scores,mx,dv);
    dt = length(film(1,1,:));
    dx = length(film(:,1,1));
    dy = length(film(1,:,1));
    %Get the time per frame
    load(string(file_path{cur})+'Time_Reference.mat');
    load(string(file_path{cur})+'Time_Groups.mat');
    t = time_ref.Y;
    t=t-t(1);
    f(cur) = uifigure('Units','Normalized','Name',string(EXP(cur)));
    tabgp(cur) = uitabgroup(f(cur),'Units','Normalized','Position',[0 0 1 1 ]);
    %Display information on each dv
    for i=1:dv
        tabs(cur,i)=uitab(tabgp(cur),'Units','Normalized','Title',"EV "+num2str(i));
        ax = axes(tabs(i),'Position',[0.01,0.95-0.4*16/9,0.4*dy/dx,0.4*16/9]);
        ax.XLim =[0 dy];
        ax.YLim = [0 dx];
        regg=col2im(squeeze(eigen(i,:)),[1 1],[dx dy]);
        imagesc(ax,regg);
        colorbar(ax);
        %Display temporal component
        ax=axes(tabs(cur,i),'Position',[0.03+0.4*dy/dx,0.7,0.44,0.2]);
        plot(ax,t,scores(i,:));
        %Determine the relevant ax limits
        avg = nanmean(scores(i,:)); 
        sig = nanstd(scores(i,:));
        miny =avg-2*sig;
        maxy = avg+2*sig;
        ax.YLim=[miny maxy];
        ax.XLim(2)=t(end);
        contrib = zeros(6,1); %to compute contribution of each phase (for pie chart)
        labels = cell(6,1);
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
        ax=axes(tabs(i),'Position',[0.05+0.4*dy/dx,0,0.4,0.6]);
        p=pie(ax,contrib); %Pie chart
        for iii=1:4
            p(2*iii - 1).FaceColor=colors(iii);
        end
        %Create buttons for Phase selection
        btn_r(i) = uibutton(tabs(cur,i),"State","Text","REM","Position",[10 30 80 20]);
        btn_nr(i) = uibutton(tabs(cur,i),"State","Text","NREM","Position",[10 5 80 20]);
        btn_w(i) = uibutton(tabs(cur,i),"State","Text","Wake","Position",[95 30 80 20]);
        btn_noise(i) = uibutton(tabs(cur,i),"State","Text","Noise","Position",[95 5 80 20]);
    end
    tabs(cur,dv+1) = uitab(tabgp(cur),'Title',string(EXP(cur)));
    %Options for algorithm to apply
    uibutton(tabs(cur,dv+1),"Text","Reconstruction","Position",[10 10 100 100],"ButtonPushedFcn",@(btn,event) pca_reconstr(eigen,scores,btn_r,btn_nr,btn_w,btn_noise,mx,phases,file_path{cur}));
    uibutton(tabs(cur,dv+1),"Text","Dimensionality","Position",[120 10 100 100],"ButtonPushedFcn",@(btn,event) pca_dim(explained,btn_r,btn_nr,btn_w,btn_noise));
    exp(cur,:) = explained(1:20);
    if bnrem
        uibutton(tabs(cur,dv+1),"Text","RIPPLES","Position",[230 10 100 100],"ButtonPushedFcn",@(btn,event) pca_ripple(btn_nr,eigen,file_path{cur},dx,dy,mx,scores,EXP));
    end
end
exp_tot = nanmean(exp,1);
exp_std = nanstd(exp,0,1);
bar(exp_tot);
hold on
er =errorbar((1:dv),exp_tot,exp_std,exp_std);
er.Color = [.8 0 .1];                            
er.LineStyle = 'none';
hold off
