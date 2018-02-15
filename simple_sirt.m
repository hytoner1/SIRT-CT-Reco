function [x, varargout] = simple_sirt(A, b, opt)
% Input: sparse system matrix A, data b.
% Output: SIRT reconstruction x.
% x = zeros(d * d, 1);

sigma = 1;
b_blur = imgaussfilt(b, sigma);

if ~isfield(opt, 'maxstep')
    opt.maxstep = 100;
end
if ~isfield(opt, 'plotFlag')
    opt.plotFlag = false;
end
if opt.plotFlag
    figure;
end
if ~isfield(opt, 'plotConv')
    opt.maxstep = false;
end
if ~isfield(opt, 'convThrs') || opt.convThrs == false
    opt.convThrs = -1e-16;
end


[Nrows, Ncols] = size(A);
x = zeros(Ncols, 1);

C = sparse(1 : Ncols, 1 : Ncols, 1 ./ sum(A));
R = sparse(1 : Nrows, 1 : Nrows, 1 ./ sum(A'));
CATR = C * A' * R;

delta = zeros(opt.maxstep,1);

for i = 1 : opt.maxstep
%     x_new = x + CATR * (b - A * x);
    
    x_blur = imgaussfilt(reshape(A*x, size(b)), sigma);
    x_new = x + CATR * (b_blur(:) - x_blur(:));
    x_new(x_new > 1) = 1; x_new(x_new < 0) = 0;
    delta(i) = sqrt(sum(x_new.^2 - x.^2));
    x = x_new;
    
    if opt.plotFlag
        imagesc(reshape(x,128,128))
        title(sprintf('Iteration %d',i));
        drawnow
    end
    if delta(i) < opt.convThrs
        break;
    end
end

if opt.plotConv
    figure; 
    plot(delta);
end

x = reshape(x, 128,128);

if max(nargout,1) == 2
    varargout = {delta};
end

end