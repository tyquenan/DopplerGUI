function [eigen,scores,mx]= pca_t(film,dv)
    dt = length(film(1,1,:));
    dx = length(film(:,1,1));
    dy = length(film(1,:,1));
    X = zeros(dx*dy,dt);
    for i =1:dt
        X(:,i)=im2col(squeeze(film(:,:,i)), [dx dy],'Distinct');
    end
    X=X';
    %% Evaluate covariance
    %x = im2col(image,[8 8],'distinct'); % convert to 8x8 blocks, each block in column
    mx = nansum(X,2)/(dx*dy);
    [eigen,scores] = pca(X);   % PCA of X returns PC and 
    eigen= eigen(:,1:dv);
    scores = scores(:,1:dv);
end
