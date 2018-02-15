%% Filesize demo for Ex 3.

%% Create a demo image

sizes = {64, 128, 256, 512, 1024, 2048};
Ps =  cellfun(@(x) phantom(x), sizes, 'UniformOutput', false);

figure(1); clf;
    for i=1:length(Ps)
        subplot( ceil(numel(sizes)/3), 3, i);
        imshow(Ps{i}, []);
        axis off;
        title(sizes{i})
    end
    
    
    %% Check the file sizes for each image, if they were isotropic 3D images

fSizes  = zeros(size(Ps));
varInfo = cell(size(Ps));

for i=1:length(Ps)
    tmp = Ps{i};
    varInfo{i} = whos('tmp');
    
    fSizes(i) = varInfo{i}.bytes * varInfo{i}.size(1);
end

%% Convert to 8-bit unsigned integers, visualize and calculate the filesizes

Ps_uint8 = cellfun(@(x) uint8(x .* 255), Ps, 'UniformOutput', false);

figure(2); clf;
    for i=1:length(Ps_uint8)
        subplot( ceil(numel(sizes)/3), 3, i );
        imshow(Ps_uint8{i}, []);
        axis off;
        title(sizes{i})
        colormap('gray')
    end
    
    %% Same filesize calculation than before
    
fSizes_uint8  = zeros(size(Ps_uint8));
varInfo_uint8 = cell(size(Ps_uint8));

for i=1:length(Ps_uint8)
    tmp = Ps_uint8{i};
    varInfo_uint8{i} = whos('tmp');
    
    fSizes_uint8(i) = varInfo_uint8{i}.bytes * varInfo_uint8{i}.size(1);
end


%% Calculate radon transform for each 

radon_times         = zeros(size(Ps));
radon_times_uint8	= zeros(size(Ps_uint8));

iradon_times        = zeros(size(Ps));
iradon_times_uint8	= zeros(size(Ps_uint8));


for i=1:length(Ps)
    tic;
    tmp = radon(Ps{i}, 1:180);
    radon_times(i) = toc;
    
    tic;
    iradon(tmp, 1:180);
    iradon_times(i) = toc;
    
    
    tic;
    tmp = radon(Ps_uint8{i}, 1:180);
    radon_times_uint8(i) = toc;
    
    tic;
    iradon(tmp, 1:180);
    iradon_times_uint8(i) = toc;
end
    
   
    %% Visualize
    
figure(3); clf; hold on;
    plot(cell2mat(sizes), radon_times .* cell2mat(sizes) ./ 60, 'b-x'); 
    plot(cell2mat(sizes), radon_times_uint8 .* cell2mat(sizes) ./ 60, 'b--x')
    
    plot(cell2mat(sizes), iradon_times .* cell2mat(sizes) ./ 60, 'r-x'); 
    plot(cell2mat(sizes), iradon_times_uint8 .* cell2mat(sizes) ./ 60, 'r--x')
    
	legend('radon double', 'radon uint8', ...
        'iradon double', 'iradon uint8',...Intersectional feminism an
        'Location', 'NW');
    ylabel('processing time (min)');
    xlabel('image side length (px)');
    title('Processing time comparison');
    
    
    
    
    
    