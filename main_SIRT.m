%% CT algebraic reco. using Simultaneous Iterative Reconstruction Technique
% Each iteration is of form 
%   v(k+1) = v(k) + CW'R (p - W'v(k)), where
%   v is reconstruction on k^th iteration, 
%   p is the original FWD projection,
%   W is the projection matrix, and
%   C and R are weight matrices (see below)

    %% Create a measurement to reconstruct

	% We want to keep axis visible
iptsetpref('ImshowAxesVisible','on')
    % Image size (square)
imS = 100;
    % Image itself
im = ones(imS);
%     im(25:75, 25:75) = 1;
    % Radon angles
theta = 0:180;
    % Radon FWD projection
[p, p_pos] = radon(im,theta);

    % Visualize
figure(1); clf;
imshow(p,[],'Xdata',theta,'Ydata',p_pos,'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gca,hot), colorbar;

    %% Reconstruct using SIRT
    
    % Start with a tabula rasa 'guess'
v = zeros(imS);
    % Crete FWD projection p_ as approximation of p
p_ = radon(v, theta);
    % Compute the difference between the two
p_diff = p - p_;

    % Introduce weight matrices 





    %% TEST_1
% Let's find the mapping matrix MAP by rotating an 'index matrix', and then checking
%  the columns

    % Create an index matrix with appropriate padding
indMat = reshape(1:imS^2, imS, imS);
    indMat = padarray(indMat, [23, 23], 'pre');
    indMat = padarray(indMat, [22, 22], 'post');

    % Allocate a map between pixel indices and sinogram values
map = cell(size(indMat, 1), length(theta));
    
	% Rotate and collect pixels affecting each 'ray'
for i = 1:length(theta)
    indMat_ = imrotate(indMat, theta(i), 'nearest', 'crop');
    
        % Collect values in each column to map
    for j = 1:size(indMat_, 1)
        tmp = indMat_(:,j);
            tmp = tmp(tmp~=0);
        map{j, i} = tmp;
    end
end


    %% TEST_1_2
% Create a sinogram on the basis of map and compare (visually) to radon

sino2 = zeros(size(map));

for i = 1 : size(map,1)
    for j = 1 : size(map,2)
        sino2(i,j) = sum( im(map{i,j}) );
    end
end

figure(2); clf; 
imshow(sino2,[],'Xdata',theta,'Ydata',p_pos,'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gca,hot), colorbar;

