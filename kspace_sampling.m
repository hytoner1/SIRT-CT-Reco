
saveFlag = 0;

im = zeros(256);

wave = sin(linspace(0,50*pi, size(im,1))) + 1;

im(:) = repmat(wave.*255, 256,1);


figure(1); clf;
    subplot(2,3,1)
    imagesc(im);
    title({'Original image', '(256x256)'})
    axis off;
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo1_a', '-dpng');
    
%%

[R, xp] = radon(im, 1:10:180);
bpR  = iradon(R, 10, 'linear','none');
bpR1 = iradon(R(:,1), 1, 'linear','none');
fbpR = iradon(R, 10);

    subplot(2,2,3);
    imagesc(R);
    title('Sinogram')
    xlabel('Angle (deg)')
    ylabel('Ray #')
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo1_b', '-dpng');
    
    subplot(2,3,2);
    imagesc(bpR1);
    axis off;
    title({'Reco from', 'First row'})
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo1_c', '-dpng');
    
    subplot(2,2,4);
    imagesc(bpR);
    title({'Reco from', 'full data'})
    axis off;
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo1_d', '-dpng');
    
    subplot(2,3,3);
    imagesc(fbpR);
    title('Filtered Reco');
    axis off;
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo1_e', '-dpng');

%%

   
freqs=linspace(-1, 1, length(xp)).';
myFilter = abs( freqs );
myFilter = repmat(myFilter, [1 size(R,2)]);

% do my own FT domain filtering
ft_R = fftshift(fft(R,[],1),1);
filteredProj = ft_R .* myFilter;
filteredProj = ifftshift(filteredProj,1);
ift_R = real(ifft(filteredProj,[],1));

figure(2); clf; 
    gs = 0.25;
    CreateAxes(2,3,1,gs);
    imagesc(im);
    axis off; axis tight
    title('Orig. image')
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_a', '-dpng');

    CreateAxes(2,3,4,0.35, [0.02,-0.1]);
    imagesc(R);
    title('Sinogram')
    xlabel('Angle (deg)')
    ylabel('Ray #')
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_b', '-dpng');
    
    CreateAxes(2,3,2,0.35);
    imagesc(log(abs(real(ft_R))));
    title('FFT of Sinogram')
    xlabel('Angle (deg)')
    ylabel('Freq');
    yticks([1, xp(end), 2*xp(end)]);
    yticklabels({xp(1), 0, xp(end)});
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_c', '-dpng');
    
    CreateAxes(2,3,5,0.35);
    plot(xp, abs(xp)./128, 'r');
    title('Ramp filter')
    xlabel('Freq')
    ylabel('Weight');
    axis tight
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_d', '-dpng');
    
    CreateAxes(2,3,3,0.35);
    imagesc(real(log(real(filteredProj))));
    title('Filtered FFT')
    xlabel('Angle (deg)')
    ylabel('Ray #')
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_e', '-dpng');
    
    CreateAxes(2,3,6,gs);
    imagesc(ift_R);
    axis off; axis tight
    title('Filtered Sino')
    
SaveCurrentFig(saveFlag, 1, './Pics', 'demo2_f', '-dpng');
    
%%

fftR = fftshift(fft(R,[],1));

figure(3); clf;
    imagesc(real(fftR));

    
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
    imagesc(R)
    subplot(2,2,2)
    imagesc(ift_R)
    subplot(2,2,3)
    imagesc(real(ft_R))
    subplot(2,2,4)
    imagesc(real(filteredProj))

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

% spect = reshape(full(LFs2) * P(:), 185,180);
spect = radon(P, 1:180);

spect = spect + (rand(size(spect))-0.5) .* 2;


figure;
    imshow(spect,[]);

%% Simple sirt

opt.maxstep = 1000;
opt.plotFlag= false;
opt.plotConv= true;
opt.convThrs= false;

tic
[reco_sirt1, delta1] = simple_sirt(LFs2, spect, opt);
toc

