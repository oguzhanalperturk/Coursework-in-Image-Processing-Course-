image1 = imread("Image1.png");
%reading image
figure;
%creating figure
subplot(1,2,1),imshow(image1);
%for showing image and its histogram on same figure
title("image1.png");
%title of image1

subplot(1,2,2),imhist(image1);
%creating histogram of image1
title("image1.png Histogram");
%title of image1 histogram


Image1Output = histeq(image1);
% I use histogram equalization to enhance the 
% image's contrast. It stretches out the intensity 
% range of the image and makes the password visible.


imwrite(Image1Output,"Image1Output.png");
% saving Image1Output.png

figure;
%creating second figure
subplot(1,2,1),imshow(Image1Output);
%showing Image1Output.png
title("Image1Output.png");
%giving title to Image1Output.png

subplot(1,2,2),imhist(Image1Output);
%showing Image1Output.png histogram
title("Image1Output.png Histogram");
%giving title to histogram




