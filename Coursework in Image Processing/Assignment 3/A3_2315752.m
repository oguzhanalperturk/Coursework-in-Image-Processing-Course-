
  image = imread("OrangesTestExample2.jpg");
  segmented_image = SegmentImage(image);

  [small_orange_count, big_orange_count] = CountOranges(segmented_image);



function [ sb ] = SegmentImage(image)

    % making the image grayscale for applying segmentation
    grayscale_image = rgb2gray(image);
    
    % Plot the histogram of an image f to see what should be threshold values
    %figure; 
    %imhist(grayscale_image);

    % Use Otsu method (graythresh built-in function) to get a threshold 
    % value between 0 and 1
        thresh = graythresh(grayscale_image);
    
    % Convert the otsu threshold value to be between 0-255
        thresh = thresh*255;
    
    % Apply bi-level thresholding (slide 35) and create new image sb

        sb=grayscale_image;
    
        sb(sb<=thresh)=0;
        sb(sb>thresh)=255;
    
    
    % Plot the results
        figure; 
        subplot(1,6,1); imshow(image);title('Original');
        subplot(1,6,2); imshow(grayscale_image);title('Grayscale');
        subplot(1,6,3); imshow(sb);title('Bilevel Thresholded');

end


function [smallObjectCount, bigObjectCount ] = CountOranges(image)

    % Image segmented using using bi-level thresholding-based segmentation. However, 
    % there are some small noises in the image. 
    % Hence, next I will make use Morphological operations and remove those
    % small noise to have a better segmentation result.
    %--------------------------------------------------------------------------

    % First I am taking iverse of the image to make imfill operation at the end.
    image = imcomplement(image);

    % Remove small pixel noises in segmented image BI using imerode built-in 
    % function and strel built-in function with 'disk' size 1. And obtain 
    % image E
        S = strel('disk',1);    % S is the structuring element
        E = imerode(image,S);   % E is the eroded version of image with S 
        
    % Plot the noise removed image 

        subplot(1,6,4);imshow(E);title('Eroded Image'); 

    % Morphology 
    % As you can see from the above figure,Erosion removes some parts in coins
    % so fill those parts of image E using imfill built-in function with option
    % 'holes' and obtain image F   

        F = imfill(E,'holes');
    
    % Plot the results to see the noise removed and filled image 

        subplot(1,6,5);imshow(F);title('Filled Image');
    
    % bwlabel returns objectCount, the number of connected objects found in the image.
    [L,objectCount] = bwlabel(F);
    
    % objectCount -> count of the objects in the image 
    % nnz(F) -> total area based on number of white pixels

    % If I divide all white pixel area to the count of the objects, the
    % number I found would be definitely between the size of the objects : 
    % small object size < P < Large object size
    % So, I can determine my P value for ismember function in that way. 

    P = nnz(F) / objectCount ;
    
    % Determine the connected components:
    CC = bwconncomp(F);

    % Compute the area of each component:
    S = regionprops(CC, 'Area');

    % Remove small objects:
    Lmat = labelmatrix(CC);
    BigObjectImage = ismember(Lmat, find([S.Area] >= P));

    % Plot the image that includes only big objects
    subplot(1,6,6);imshow(BigObjectImage);title('Big Objects Only');

    % Counting the remaining big objects
    [L,bigObjectCount] = bwlabel(BigObjectImage);

    % Small objects = All objects - Big Objects
    smallObjectCount = objectCount - bigObjectCount;

    figure;
    subplot(1,3,1);imshow(F);title("Segmented Image");
    subplot(1,3,2);text( 0.4, 0.5, string(bigObjectCount), 'fontsize', 48 );title("Big Object Count");
    subplot(1,3,3);text( 0.4, 0.5, string(smallObjectCount), 'fontsize', 48 );title("Small Object Count");

end







