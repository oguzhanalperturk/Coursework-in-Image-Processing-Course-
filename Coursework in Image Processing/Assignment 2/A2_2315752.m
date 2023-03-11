% Oğuzhan Alpertürk
% 2315752

clear all;clc;

%% Question1
% read the image
I1 = imread("noisy1.png");

% I investigate the image by looking at its histogram.
% figure;
% imhist(I1)
% I also took a block and investigate its histogram as well
% figure;
% imshow(I1(1:100,1:100));
% imhist(I1(1:100,1:100));
% After looking at these histograms and comparing the image with the images 
% in the lecture slides, the noise in the image is not an additive noise type.

% take the fourier transform of the image to investigate it in the
% frequency domain. 

% Compute the shifted DFT of the image using functions fft2 and fftshift
F= fftshift(fft2(I1));

% figure;
% imshow(log(abs(F)),[]);
% after looking at it in frequency domain, I decided that its a periodic
% noise.

% After analyzing the noise in the frequency domain, I decided that
% gaussian lowpass filter can be used for handling the noise in the image.

% Create the gaussian lowpass filter.
D0=100; n=2;
[M, N] = size(I1);
[Hl,Hh] = LowAndHighPassFilters('gaussian',M,N,D0,n);

% Filter the image by multiplying the filter with shifted DFT of the image
G = Hl .* F;

% Compute the inverse DFT using ifft2 and abs functions
G2 = abs(ifft2(G));

% Convert double to image using uint8 function 
G3= uint8(255 * mat2gray(G2));

% Creating/Writing recovered image
imwrite(G3, 'recovered1.png');

% find edges of original image
Edge_original = findEdges(I1);

% find edges of denoised image 
Edge_denoised = findEdges(G3);

% subtract them
edge_difference = Edge_original - Edge_denoised;

figure;
imshow(edge_difference);


%% Question2
% read the image
I2 = imread("noisy2.png");

% I investigate the image by looking at its histogram.
% figure;
% imhist(I2)
% I also took a block and investigate its histogram as well
% figure;
% imshow(I2(1:100,1:50));
% imhist(I2(1:100,1:50));
% After looking at these histograms and comparing the image with the images 
% in the lecture slides, the noise in the image is not an additive noise type.

% Compute the shifted DFT of the image using functions fft2 and fftshift
F= fftshift(fft2(I2));

% figure;
% imshow(log(abs(F)),[]);
% after looking at it in frequency domain, I decided that its a periodic
% noise.

% After analyzing the noise in the frequency domain, I decided that
% ideal bandpass filter can be used for handling the noise in the image.
 
% Create the ideal bandpass filter.
% My aim is to cover noises in the frequency domain and I chose the needed
% values in that way.

D0=60; n=2;
[M, N] = size(I2);
[Hbp,Hbr] = BandPassAndRejectFilters('ideal',M,N,D0,n,15);

% Filter the image by multiplying the filter with shifted DFT of the image
G =  Hbr .* F;

% Compute the inverse DFT using ifft2 and abs functions
G2 = abs(ifft2(G));

% Convert double to image using uint8 function 
G3= uint8(255 * mat2gray(G2));

% Creating/Writing recovered image
imwrite(G3, 'recovered2.png');

% find edges of original image
Edge_original = findEdges(I2);

% find edges of denoised image 
Edge_denoised = findEdges(G3);

% subtract them
edge_difference = Edge_original - Edge_denoised;

figure;
imshow(edge_difference);


%% Question3
% read the image
I3 = imread("noisy3.tif");

% I investigate the image by looking at its histogram.
% figure;
% imhist(I3)
% I also took a block and investigate its histogram as well
% figure;
% imshow(I3(1:100,1:50));
% imhist(I3(1:100,1:100));
% After looking at these histograms and comparing the image with the images 
% in the lecture slides, the noise in the image is not an additive noise type.

% Compute the shifted DFT of the image using functions fft2 and fftshift
F= fftshift(fft2(I3));

% figure;
% imshow(log(abs(F)),[]);
% after looking at it in frequency domain, I decided that its a periodic
% noise.

% After analyzing the noise in the frequency domain, I decided that
% combination of the ideal notchreject filter can be used for handling 
% the noise in the image.

% Create the ideal notchreject filter.
% My aim is to cover noises in the frequency domain and I chose the needed
% values in that way.
% Because of the noise type (8 spot should be masked) , I combined 4 notchreject filter.
D0=12; n=2;
[M, N, dim] = size(I3);
[Hnp1,Hnr1] = NotchPassAndRejectFilters('ideal',M,N,D0,40,32,n);
[Hnp2,Hnr2] = NotchPassAndRejectFilters('ideal',M,N,D0,-40,32,n);
[Hnp3,Hnr3] = NotchPassAndRejectFilters('ideal',M,N,D0,80,32,n);
[Hnp4,Hnr4] = NotchPassAndRejectFilters('ideal',M,N,D0,-80,32,n);

filter = Hnr1 .* Hnr2 .* Hnr3 .* Hnr4;

% Filter the image by multiplying the filter with shifted DFT of the image
G =  filter .* F;

% Compute the inverse DFT using ifft2 and abs functions
G2 = abs(ifft2(G));

% Convert double to image using uint8 function
G3 = uint8(255 * mat2gray(G2));

% Creating/Writing recovered image
imwrite(G3, 'recovered3.png');

% find edges of original image
Edge_original = findEdges(I3);

% find edges of denoised image 
Edge_denoised = findEdges(G3);

% subtract them
edge_difference = Edge_original - Edge_denoised;

figure;
imshow(edge_difference);


function [Hl,Hh] = LowAndHighPassFilters(TYPE,M,N,D0,n)
%Computes frequency domain lowpass and highpass filters
%   H = LowAndHighPassFilters(TYPE, M, N, D0, n) creates a lowpass, Hl, and
%   a highpass filter, Hh, of the specified TYPE and size (M-by-N). 
%
%   Valid values for TYPE, D0, and n are:
%   'ideal'       Ideal filter with cutoff frequency D0.
%
%   'butterworth' Butterworth filter of order n, and cutoff D0.
%
%   'gaussian'    Gaussian filter with cutoff frequency DO.

% Compute the distances D(U, V)
    for u=1:M
        for v=1:N
         D(u,v)=((u-(M/2))^2 + (v-(N/2))^2 )^(1/2);
         end
    end
    
% Compute the filter   
    switch TYPE
        case 'ideal' 
             Hl=zeros(M,N);
             Hl(D<=D0)=1;
             Hh=1-Hl;
        case 'butterworth'
             Hl = 1./(1 + (D./D0).^(2*n)); 
             Hh=1-Hl;
        case 'gaussian'
             Hl =exp(-(D.^2)./(2*(D0^2)));
             Hh=1-Hl;
    end
end


function [Hbp,Hbr] = BandPassAndRejectFilters(TYPE,M,N,D0,n,W)
%Computes frequency domain bandpass and bandreject filters.
%   H = BandPassAndRejectFilters(TYPE, M, N, D0, n) creates a bandpass, Hbp, 
%   and a bandreject filter, Hbr, of the specified TYPE and size (M-by-N). 
%
%   Valid values for TYPE, D0, and n are:
%   'ideal'         Ideal filter with cutoff frequency D0 and width W.
%
%   'butterworth'   Butterworth filter of order n, cutoff D0 and width W.
%
%   'gaussian'      Gaussian filter with cutoff frequency D0 and width W.

% Compute the distances D(U, V)
    for u=1:M
        for v=1:N
         D(u,v)=((u-(M/2))^2 + (v-(N/2))^2 )^(1/2);
         end
    end
    
% Compute the filter   
    switch TYPE
        case 'ideal' 
            % bandpass Hbp and a bandreject Hbr filter 
            Hbr = ones(M,N);
            Hbr(D >= (D0 - (W/2)) & D <= D0 + (W/2)) = 0;
            Hbp = 1 - Hbr;
        case 'butterworth'
            % bandpass Hbp and a bandreject Hbr filter
            Hbr = 1./(1 + (((D.*W) ./ (D.^2 - D0.^2)).^(2*n)));
            Hbp = 1 - Hbr;
        case 'gaussian'
            % bandpass Hbp and a bandreject Hbr filter 
            Hbr = 1 - exp(-((D.^2 - D0.^2) ./ (D.*W)).^2);
            Hbp = 1 - Hbr;
    end
end


function [Hnp,Hnr] = NotchPassAndRejectFilters(TYPE,M,N,D0,uk,vk,n)
%Computes frequency domain notchpass and notchreject filters.
%   H = NotchPassAndRejectFilters(TYPE, M, N, D0, n) creates a notchpass,  
%   Hnp,and a notchreject filter, Hnr, of the specified TYPE and 
%   size (M-by-N). 
%
%   Valid values for TYPE, D0, and n are:
%   'ideal'         Ideal  filter with cutoff frequency D0 and centers(uk,vk).
%
%   'butterworth'   Butterworth filter of order n, cutoff D0 and centers(uk,vk).
%
%   'gaussian'      Gaussian filter with cutoff frequency DO and centers(uk,vk).

% Define D(u,v)
    for u=1:M
        for v=1:N
         Dkp(u,v)=((u-(M/2)-uk)^2 + (v-(N/2)-vk)^2 )^(1/2);
         Dkn(u,v)=((u-(M/2)+uk)^2 + (v-(N/2)+vk)^2 )^(1/2);
         end
    end
    
% Compute the filter   
    switch TYPE
        case 'ideal' 
            % notchpass Hnp and a notchreject Hnr filter
            Hnr = ones(M,N);
            Hnr(Dkp <= D0 | Dkn <= D0) = 0;
            Hnp = 1-Hnr;
        case 'butterworth'
            % notchpass Hnp and a notchreject Hnr filter
            Hnr = 1-((1./(1 + (Dkp./D0).^(2*n))) .* 1./(1 + (Dkn./D0).^(2*n)));
            Hnp = 1-Hnr;
        case 'gaussian'
            % notchpass Hnp and a notchreject Hnr filter
            Hnr = (1-exp(-(Dkp.^2)./(2*(D0^2)))) .* (1-exp(-(Dkn.^2)./(2*(D0^2))));
            Hnp = 1-Hnr;
    end
end



function [G4] = findEdges(I)

% I used sobel filters for finding the edges

sobel_x = [-1,-2,-1;
            0,0,0;
            1,2,1];

sobel_y = [-1,0,1;
           -2,0,2;
           -1,0,1];
[M,N] = size(I);

 Hx = (fftshift(abs(fft2(sobel_x,M,N)))); 
 Hy = (fftshift(abs(fft2(sobel_y,M,N))));

% Compute the shifted DFT of the image using functions fft2 and fftshift
 F= fftshift(fft2(I));

% Filter the image with multiplying the filter with shifted DFT of the image
 G_x = F .* Hx;
 G_y = F .* Hy;


% Compute the inverse DFT using ifft2 and abs functions
 G1 = uint8(abs(ifft2(G_x))); 
 G2 = uint8(abs(ifft2(G_y))); 

 G3 = G1 + G2;

% Convert double to image using uint8 function    
 G4 = uint8(255 * mat2gray(G3));

end


