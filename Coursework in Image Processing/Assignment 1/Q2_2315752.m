
image = imread("Image2.png");
%read the image

figure;
subplot(1,4,1),imshow(image);
title("Image2.png");
%showing original image 

% noise type : salt and pepper 
% I used median filter to denoise the image. 
denoised_im = applyMedianFilter(image);

imwrite(denoised_im, "Image2Output.png");
%save the denoised image
subplot(1,4,2),imshow(denoised_im);
%show denoised image
title("Image2Output.png");


edge_image = findEdges(image);
%find edges of original image with sobel filters
subplot(1,4,3),imshow(edge_image);
%show original edge image
title("Edges of Image2.png");


edge_denoised_im = findEdges(denoised_im);
% find edges of denoised image
subplot(1,4,4),imshow(edge_denoised_im);
% show denoised edge image
title("Edges of Image2Output.png");
% The edges were preserved.




function [ conv_image ] = calculateConv_w3x3filter(image,filter)

[height, width] = size(image);
conv_image = zeros(height,width);

% rotate filter by 180 degrees (convolution)
filter = [filter(1,1),filter(2,1),filter(3,1);
          filter(1,2), filter(2,2),filter(3,2);
          filter(1,3),filter(2,3),filter(3,3)];

% reaching all pixels with for loop. 
% multiply values which overlap the image with the corresponding image intensity value.
% Then sum them up and write it into the center pixel.

    for i=1:height
        for j=1:width
            if(i==1 && j == 1)
                conv_image(i,j) = (image(i,j)*filter(2,2))+image(i,j+1)*filter(2,3)+image(i+1,j)*filter(3,2) + image(i+1,j+1)*filter(3,3);  
            
            elseif(i == 1 && j == width)
                conv_image(i,j) = image(i,j-1)*filter(2,1) + (filter(2,2)*image(i,j)) + image(i+1,j)*filter(3,2) + image(i+1,j-1)*filter(3,1);
          
            elseif(i == 1)
                conv_image(i,j) = image(i,j-1)*filter(2,1) + (filter(2,2)*image(i,j)) + image(i,j+1)*filter(2,3) + image(i+1,j-1)*filter(3,1) + image(i+1,j)*filter(3,2) + image(i+1,j+1)*filter(3,3);
            
            elseif(i==height && j ==1)
                conv_image(i,j) = image(i-1,j)*filter(1,2) + (image(i,j)*filter(2,2)) + image(i,j+1)*filter(2,3) + image(i-1,j+1)*filter(1,3);
    
            elseif(i==height && j ==width)
                conv_image(i,j) = image(i-1,j)*filter(1,2) + image(i,j-1)*filter(2,1) + (filter(2,2)*image(i,j)) + image(i-1,j-1)*filter(1,1);
            
            elseif(i==height)
                conv_image(i,j) = image(i-1,j-1)*filter(1,1) + image(i-1,j)*filter(1,2) + image(i-1,j+1)*filter(1,3) + image(i,j-1)*filter(2,1) +  image(i,j)*filter(2,2) + image(i,j+1)*filter(2,3);
            
            elseif(j==width)
                conv_image(i,j) = image(i-1,j-1)*filter(1,1) + image(i,j-1)*filter(2,1) + image(i+1,j-1)*filter(3,1) + image(i-1,j)*filter(1,2) + image(i,j)*filter(2,2) + image(i+1,j)*filter(3,2);
            
            elseif(j == 1)
                conv_image(i,j) = image(i-1,j)*filter(1,2) + image(i,j)*filter(2,2) + image(i+1,j)*filter(3,2) + image(i-1,j+1)*filter(1,3) + image(i,j+1)*filter(2,3) + image(i+1,j+1)*filter(3,3);
            else
                conv_image(i,j) = image(i-1,j-1)*filter(1,1) + image(i-1,j)*filter(1,2) + image(i-1,j+1)*filter(1,3) + image(i,j-1)*filter(2,1) + image(i,j)*filter(2,2) + image(i,j+1)*filter(2,3) + image(i+1,j-1)*filter(3,1) + image(i+1,j)*filter(3,2) + image(i+1,j+1)*filter(3,3);
            end
        end
    end         
end


function[ edge_im ] = findEdges(image)

% sobel filters for detecting edges
    
    sobel_x = [-1,-2,-1;
               0,0,0;
               1,2,1];

    sobel_y = [-1,0,1;
               -2,0,2;
               -1,0,1];
    
    % for visualization purposes:
    % I used abs -> taking the absolute value of the image
    % uint8 -> turn into the image from double to integer 

    % edges in x dimention
    conv_tif6_x = calculateConv_w3x3filter(double(image),sobel_x);
    %convolution of the image with sobel_x
    conv_tif6_x = abs(conv_tif6_x);
    conv_tif6_x = uint8(conv_tif6_x);
    
    %edges in y dimention
    conv_tif6_y = calculateConv_w3x3filter(double(image),sobel_y);
    %convolution of the image with sobel_y
    conv_tif6_y = abs(conv_tif6_y);
    conv_tif6_y = uint8(conv_tif6_y);
    
    % sum them up for detecting all edges
    edge_im = (conv_tif6_y + conv_tif6_x);
end


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










