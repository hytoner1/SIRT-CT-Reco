%Author: Stefano "grog" Ghio, Michele Bozzano, Bruno Mazzarello
%Released under CreativeCommons and GPLv2 licenses

%SIRT - Simultaneous Iterative Reconstruction Technique
%Formula is:
%f_(k+1) = f_k + At (g - A f_K)
%f_k is solution at k-th iteration, at first iteration it is our guess
%g is image sinogram
%A f_k is Radon transform of f_k
%At is inverse Radon of its argument

%This example shows how to scan an image with a limited angle, obtain the sinogram, add noise to make things worse, and reconstruct the initial image with great accuracy
%You may want to remove/edit some sections if you intend to apply this to your use case:
%- skip/correct normalization
%- skip noise
%- use FBP for the iradon

%clear matlab environment
clear all;
close all;
theta = [0:179]; %projection angle 180 degrees
F = phantom(128); %create new test phantom 128x128 pixels aka alien. You can use your own image here if you want. It was not tested on non-grayscale images
figure, imshow(F,[]),title('Alien vs Predator'); %show alien
S1 = sum(sum(F));%calculate pixels sum on F
R = radon(F,theta);%apply Radon aka scan the alien to create a sinogram
figure, imshow(R,[]),title('Sinoalien');%show sinogram
%add image noise to the sinogram. You can skip this part, it's just to showcase the algorithm, it's not obviously needed
maximum=max(max(R));
minimum=min(min(R));
R=(R-minimum)/(maximum-minimum);
%normalize between 0 and 1. You can tune this as per your needs
fact=1001;%set the number of X-rays, the higher the better (and deadlier)
R=(fact/10^12)*R;
Noised=imnoise(R,'poisson');%add Poissonian noise
figure, imshow(Noised,[]),title('Dirty alien');%show noisy sinogram
%This part below you actually need it! Edit accordingly to your case
%At applied to g aka noisy alien
%Check the parameters here if you edited the code above! You might also want to tune the parameters for the iradon. Here we're showing that this works even when the FBP is not available!
At = iradon(Noised,theta,'linear', 'none', 1,128); %reconstruct noisy alien.
figure, imshow(At,[]),title('At G'); %show noisy alien
%algorithm starts crunching here
S2 = sum(sum(At)); %calculate pixels sum on At
At = (At/S2)*S1; %normalize At so that pixel counts match. Edit as per your needs
n = 100;%iterations. Might be more, there's always a limit over which it doesn't make sense to keep iterating though!
Fk = At;%Matrix Fk is our solution at the k-th step, now it is our initial guess
for  k=1:n
    t = iradon(radon(Fk,theta),theta, 'linear', 'none', 1,128);% reconstruct alien using Fk unfiltered sinogram. Maybe use FBP if you want here
    %normalize. Again, as per your needs
    St = sum(sum(t));
    t = (t/St)*S1;
    %update using (At g - At A f_k) 
    %new Fk = Fk + difference between reconstructed aliens At_starting - t_previuous_step
    Fk = Fk + At - t;
	%remember that our alien is normalized between 0 and 1. Might not be your case!
    %delete values <0 aka not real! Might not be your case!
    index = Fk<0;
    Fk(index)=0;
    %delete values >1 aka not real! Might not be your case!
    index = find(Fk>1);
    Fk(index)=1;
    %show reconstruction progress every 10 steps
    if(mod(k,10)== 0)
    figure,imshow(Fk,[]),title('Fk');
    end
    %calculate step improvement between steps. Tune as per your needs
    Arrerr(k) = sum(sum((F - Fk).^2));
    %stop when there is no more improvement. Tune as per your needs
    if((k>2) &&(Arrerr(k)>Arrerr(k-1)))
       break;
    end
end
figure, plot(Arrerr),title('Arrerr');%show error improvement over all iterations
