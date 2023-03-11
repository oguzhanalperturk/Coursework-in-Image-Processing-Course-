
image = imread("image3.png");
% read the image3.png
figure;
subplot(1,2,1),imshow(image);
title("Image3.png");
% show image3.png

% treshold the image
% I chose the value of 142 because it makes median filter to denoise all the
% small white dots.
tresholded_image = image>142;

% I applied median filter to tresholded image
% The reason is to clean the small white dots in the tresholded image
tresholded_image = applyMedianFilter(tresholded_image);

% After applied median filter, there is only a big white star visible.
% but the star is white (tresholded)
% I used "detectBrightestStar" function to change white intensity value with the
% original intensity value of the star.
brightest_star_im = detectBrightestStar(image, tresholded_image);

imwrite(brightest_star_im, "Image3Output.png");
% save Image3Output.png

subplot(1,2,2),imshow(brightest_star_im);
title("Image3Output.png");
% show Image2Output.png




function [denoised_im] = applyMedianFilter(image) 
    
%median_filter=[image(i-1,j-1), image(i-1,j), image(i-1,j+1);
%               image(i,j-1),   image(i,j),   image(i,j+1);
%               image(i+1,j-1), image(i+1,j), image(i+1,j+1)];

denoised_im = image;
[height, width] = size(denoised_im);
%getting height and width of the denoised image

% reaching all pixels with for loop. 
% put into the values which overlap the image in a list and then sort them.
% Then choose the middle value (median) and put it into the image pixel which overlaps the center of the filter 

    for i=1:height
        for j=1:width
            if(i==1 && j == 1)
                median_filter = sort([image(i,j), image(i,j+1), image(i+1,j), image(i+1,j+1)]);  
                denoised_im(i,j) = (median_filter(2) + median_filter(3)) / 2;
    
            elseif(i == 1 && j == width)
                median_filter = sort([image(i,j-1), image(i,j), image(i+1,j-1), image(i+1,j)]);
                denoised_im(i,j) = (median_filter(2) + median_filter(3)) / 2;
          
            elseif(i == 1)
                median_filter = sort([image(i,j-1), image(i,j), image(i,j+1), image(i+1,j-1), image(i+1,j), image(i+1,j+1)]);
                denoised_im(i,j) = (median_filter(3) + median_filter(4)) / 2;
    
            elseif(i==height && j ==1)
                median_filter = sort([image(i-1,j), image(i-1,j+1), image(i,j), image(i,j+1)]);
                denoised_im(i,j) = (median_filter(2) + median_filter(3)) / 2;
    
            elseif(i==height && j ==width)
                median_filter = sort([image(i-1,j-1), image(i-1,j), image(i,j-1), image(i,j)]);
                denoised_im(i,j) = (median_filter(2) + median_filter(3)) / 2;
    
            elseif(i==height)
                median_filter = sort([image(i-1,j-1), image(i-1,j), image(i-1,j+1), image(i,j-1), image(i,j), image(i,j+1)]);
                denoised_im(i,j) = (median_filter(3) + median_filter(4)) / 2;
    
            elseif(j==width)
                median_filter = sort([image(i-1,j-1), image(i,j-1), image(i+1,j-1), image(i-1,j), image(i,j), image(i+1,j)]);
                denoised_im(i,j) = (median_filter(3) + median_filter(4)) / 2; 
                
            elseif(j == 1)
                median_filter = sort([image(i-1,j), image(i,j), image(i+1,j), image(i-1,j+1), image(i,j+1), image(i+1,j+1)]);
                denoised_im(i,j) = (median_filter(3) + median_filter(4)) / 2;             
    
            else
                median_filter=sort([image(i-1,j-1),image(i-1,j),image(i-1,j+1),image(i,j-1),image(i,j),image(i,j+1),image(i+1,j-1),image(i+1,j),image(i+1,j+1)]);
                denoised_im(i,j) = median_filter(5);
            end
        end
    end
end



% The function changes the white values in the tresholded image with the real (original) intensity values in the original image. 

function [ brightest_star_im ] = detectBrightestStar(image, tresholded_image)

[height, width] = size(image);
% getting height and width
brightest_star_im = image;

    for i=1:height
        for j=1:width
            if(tresholded_image(i,j) == 0)  % if value is 0 then it remains the same
                brightest_star_im(i,j) = 0;
            else
                brightest_star_im(i,j) = image(i,j);  % if the value is white, then change it with the original value
            end
        end
    end         
end


