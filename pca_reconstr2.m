    function pca_reconstr2(eigen,scores,EV_rem,EV_nrem,EV_wake,EV_noiseless,mx,phases,curpath)
    load(string(curpath)+"Doppler.mat");
    film = Doppler_film;
    dt = length(film(1,1,:));
    dx = length(film(:,1,1));
    dy = length(film(1,:,1));
    %%
    films = zeros(1,dt,dx*dy);
    disp([length(EV_rem) length(EV_nrem) length(EV_wake) length(find([noise.Value] == 1))]);
    eigen_nless=eigen(EV_noiseless,:);
    scores_nless = scores(EV_noiseless,:);
    X_nless = scores_nless'*eigen_nless;
    for i=1:dx*dy
        X_nless(:,i) = X_nless(:,i)+mx;
    end
    figure("Units","normalized");
    for i=1000:dt
        ax= subplot('position',[(1-0.5*16/9)/2, 0.05,0.5*128/101,0.5*16/9]);
        imagesc(col2im(squeeze(X_nless(i,:)),[1 1],[dx dy]));
        title(sprintf('Film Noiseless Image %d/%d  phase %s',[i dt phases{i}]));
        colorbar(ax);
        ax.CLim=[-100,100];
        pause(0.05);
    end
    eigen_rem=eigen(EV_rem,:);
    scores_rem = scores(EV_rem,:);
    X_rem = scores_rem'*eigen_rem;
    for i=1:dx*dy
        X_rem(:,i) = X_rem(:,i)+mx;
    end
    figure("Units","normalized");
    for i=1:dt
        subplot('position',[(1-0.5*16/9)/2, 0.05,0.5*128/101,0.5*16/9]);
        imagesc(col2im(squeeze(X_rem(i,:)),[1 1],[dx dy]));
        title(sprintf('Film REM Image %d/%d  phase %s',[i dt phases{i}]));
        colorbar;
        pause(0.05);
    end    
    eigen_nrem=eigen(EV_nrem,:);
    scores_nrem = scores(EV_nrem,:);
    X_nrem = scores_nrem'*eigen_nrem;
    for i=1:dx*dy
        X_nrem(:,i) = X_nrem(:,i)+mx;
    end
    figure("Units","normalized");
    for i=1:dt
        subplot('position',[(1-0.5*16/9)/2, 0.05,0.5*128/101,0.5*16/9]);
        imagesc(col2im(squeeze(X_nrem(i,:)),[1 1],[dx dy]));
        title(sprintf('Film NREM Image %d/%d  phase %s',[i dt phases{i}]));
        colorbar;
        pause(0.05);
    end    
    eigen_wake=eigen(EV_wake,:);
    scores_wake = scores(EV_wake,:);
    X_wake = scores_wake'*eigen_wake;
    for i=1:dx*dy
        X_wake(:,i) = X_wake(:,i)+mx;
    end
    figure("Units","normalized");
    for i=1:dt
        subplot('position',[(1-0.5*16/9)/2, 0.05,0.5*128/101,0.5*16/9]);
        imagesc(col2im(squeeze(X_wake(i,:)),[1 1],[dx dy]));
        title(sprintf('Film wake Image %d/%d  phase %s',[i dt phases{i}]));
        colorbar;
        pause(0.05);
    end    
    quit
end