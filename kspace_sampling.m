

im = uint8(zeros(256,256));

wave = sin(linspace(0,400*pi, size(im,1))) + 1;

im(:) = repmat(wave.*255, 256,1);


% figure; imshow(im,[])

%%
im_fft = (fft2(P));

figure; plot(abs(im_fft(2,1:128)).*linspace(0,1,128))


%% 

P = phantom(128); 
    imshow(P);
    
%%

[R, xp] = radon(P, 1:180);

figure; imshow(R,[])
figure; plot(R(:,1));

r1 = R(:,40);
r1_ft = fft(r1);

figure; plot(abs(fftshift(r1_ft)))

%%

freqs=linspace(-1, 1, length(xp)).';
myFilter = abs( freqs );
myFilter = repmat(myFilter, [1 size(R,2)]);

% do my own FT domain filtering
ft_R = fftshift(fft(R,[],1),1);
filteredProj = ft_R .* myFilter;
filteredProj = ifftshift(filteredProj,1);
ift_R = real(ifft(filteredProj,[],1));

figure;
    subplot(2,2,1)
    imshow(R,[])
    subplot(2,2,2)
    imshow(ift_R,[])
    subplot(2,2,3)
    imshow(ft_R,[])
    subplot(2,2,4)
    imshow(filteredProj,[])

%%

R_reco = iradon(R, 1:180,'linear','none');
R2_reco = iradon(ift_R, 1:180,'linear','none');

figure;
    subplot(1,2,1)
    imshow(R_reco,[])
    subplot(1,2,2)
    imshow(R2_reco,[])

%%


figure(2); clf;
for i=1:180
    R_reco = iradon(R(:,1:i), 1:i,'linear','none');
    
    R2_reco = iradon(ift_R(:,1:i), 1:i,'linear','none');
    
    subplot(1,2,1);
    imshow(R_reco, []);    
    subplot(1,2,2);
    imshow(R2_reco, []);
    title(sprintf('No. angles: %d', i));
%     drawnow
pause
end

%% Lead field pööpöilyä

LFs = repmat({sparse(1)},length(P(:)),1);

tic
for j = 1:length(P(:))
	im = zeros(size(P,1), 128);
    im(j) = 1;
    tmp = radon(im);
    LFs{j} = sparse(tmp(:));
end

toc

LFs2 = sparse(cell2mat(LFs'));

%% Apply lead field to the phantom image

spect = reshape(full(LFs2) * P(:), 185,180);

figure;
    imshow(spect,[]);

%% Simple sirt

opt.maxstep = 1000;
opt.plotflag= false;
opt.plotConv= true;
opt.convThrs= false;

tic
reco_sirt = simple_sirt(LFs2, R, opt);
toc

