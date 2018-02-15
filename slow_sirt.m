function [reco, varargout] = slow_sirt(sinogram, map, imS)
%SLOW_SIRT computes angle-by-angle iterative reconstruction
%   Using the method detailed in Webb; The Physics of Medical Imaging
% INPUTS:
%   sinogram - a measurement matrix generated e.g. with radon()
%   map - mapping matrix that links sinogram pixels to patient space
% OUTPUTS:
%   reco - the reconstructed image

    % Tabula rasa, on top of which we project
reco = zeros(imS);

reco_order = randperm(size(map,2));
iter = 0;

figure(1); clf;
    imagesc(reco_order);
    title(sprintf('Reco after %d iterations', iter));
    
for i = 1:length(reco_order)
    iter = iter+1;
    for j = 1:size(map,1)
        reco(map{j,reco_order(i)}) =...
            reco(map{j,reco_order(i)}) + sinogram(j,reco_order(i));
    end

    reco_fwd = zeros(size(sinogram,1));
    for k=1:size(map,1)
        reco_fwd(k) = sum(reco(map{k, reco_order(i+1)}));
    
    imagesc(reco);
    title(sprintf('Reco after %d iterations', iter));
    drawnow;
end
    
    
end

