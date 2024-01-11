function [eigen,scores,mx] = pca_t_cov_alt(film,eigen,scores,mx,dv)
    dt = length(film(1,1,:));
    dx = length(film(:,1,1));
    dy = length(film(1,:,1));
    eigen=eigen';
    scores=scores';
end
