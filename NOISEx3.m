%%ADDING NOISE TO THE 3 SEGMENTATION METHODS
%% NOISE & MANUAL THRESHOLDING+HISTOGRAM EQUALIZATION
clear all;
close all;
clc;

basePath = 'C:\Users\micol\OneDrive\Documenti\MATLAB\MI\PROJECT\Segmentation\Patient 0\';
timeFrames = {'T_0','T_10','T_20', 'T_30', 'T_40','T_50','T_60','T_70','T_80', 'T_90'};

for t = 1:length(timeFrames)
    timeFrame = timeFrames{t};
    folderPath = fullfile(basePath, timeFrame, 'CT');
    filePattern = fullfile(folderPath, '*.dcm');
    dicomFiles = dir(filePattern);

    % ROI
    fullFileName = fullfile(folderPath, dicomFiles(1).name);
    dicomImage = dicomread(fullFileName);
    imshow(dicomImage, []);
    [croppedImage, rect] = imcrop;
    waitforbuttonpress;
    rect = round(rect);
    close all;

    % Array 3D for memorizing the volume later
    vol = zeros(rect(4), rect(3), length(dicomFiles), 'logical');

    % Processing the DICOM files
    for k = 1:length(dicomFiles)
        fullFileName = fullfile(folderPath, dicomFiles(k).name);
        fprintf(1, 'Reading %s\n', fullFileName);
        dicomImage = dicomread(fullFileName);

        croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));

        % Add Gaussian noise to the image
        noisy_image = imnoise(dicomImage,'gaussian',0,10^-4);

        croppedImage = noisy_image(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));
        croppedImageDouble = double(croppedImage);
        normalizedCroppedImage = (croppedImageDouble - min(croppedImageDouble(:))) / (max(croppedImageDouble(:)) - min(croppedImageDouble(:)));
        eqCroppedImage = histeq(normalizedCroppedImage);
        manualThreshold = 0.3
        Segmented_Lungs = imbinarize(eqCroppedImage, manualThreshold);
      
        Segmented_Lungs = ~Segmented_Lungs;
        
        Lungs = bwareafilt(Segmented_Lungs, 2);
        Lungs = bwareaopen(Lungs, 500); 
        Lungs = imfill(Lungs, 'holes');
        Lungs = bwareaopen(Lungs, 5000); 

       
        vol(:,:,k) = Lungs;
        
        % Plotting
        if k == 1
            figure;
            sgtitle('Segmentation with Manual Thresholding+ Histogram Equalization')

            subplot(2, 2, 1);
            imshow(Segmented_Lungs);
            title('Segmented Lungs');
            
            subplot(2, 2, 2);
            imshow(eqCroppedImage);
            title('Equalized Cropped Image');

            subplot(2, 2, 3);
            imshow(~Lungs); 
            title('Cleaned Lungs');

            subplot(2, 2, 4);
            imhist(eqCroppedImage);
            title('Histogram of Equalized Image');
        end
    end

    % 3D Volume
    if exist('volshow', 'file')
        figure;
        volshow(vol);
    end
end


%% NOISE & OTSU METHOD

clear all;
close all;
clc;

basePath = 'C:\Users\micol\OneDrive\Documenti\MATLAB\MI\PROJECT\Segmentation\Patient 0\';
timeFrames = {'T_0','T_10','T_20', 'T_30', 'T_40','T_50','T_60','T_70','T_80', 'T_90'};

for t = 1:length(timeFrames)
    timeFrame = timeFrames{t};
    folderPath = fullfile(basePath, timeFrame, 'CT');
    filePattern = fullfile(folderPath, '*.dcm');
    dicomFiles = dir(filePattern);

    % ROI
    fullFileName = fullfile(folderPath, dicomFiles(1).name);
    dicomImage = dicomread(fullFileName);
    imshow(dicomImage, []);
    [croppedImage, rect] = imcrop;
    waitforbuttonpress;
    rect = round(rect);
    close all; 

    % Array 3D for memorizing the volume later
    vol = zeros(rect(4), rect(3), length(dicomFiles), 'logical'); 

    % Processing the DICOM files
    for k = 1:length(dicomFiles)
        fullFileName = fullfile(folderPath, dicomFiles(k).name);
        fprintf(1, 'Reading %s\n', fullFileName);
        dicomImage = dicomread(fullFileName);

        % Add Gaussian noise to the image
        noisy_image = imnoise(dicomImage,'gaussian',0,10^-4);
        
        croppedImage = noisy_image(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));
        level = graythresh(croppedImage);
        Segmented_Lungs = imbinarize(croppedImage, level);
        Lungs = bwareafilt(~Segmented_Lungs, 2);
        Lungs = imfill(Lungs, 'holes');
        Lungs = bwareaopen(Lungs, 500);  
        Lungs = bwareaopen(Lungs, 5000); 

        % Plotting
    if k == 1
        figure;
        subplot(1, 4, 1);
        imshow(croppedImage);
        title('Cropped Image');

        subplot(1, 4, 2);
        imshow(Segmented_Lungs);
        title('Segmented Lungs');

        subplot(1, 4, 3);
        imshow(~Lungs); 
        title('Cleaned Lungs');

        subplot(1, 4, 4);
        imhist(croppedImage);
        title('Histogram of Cropped Image');
    end

        % Applying the mask to each slice
        vol(:,:,k) = Lungs;
    end

    % 3D Volume
    if exist('volshow', 'file')
        figure;
        volshow(vol);
    end
end



%% NOISE & GRADIENT METHOD

clear all;
close all;
clc;

basePath = 'C:\\Users\\micol\\OneDrive\\Documenti\\MATLAB\\MI\\PROJECT\\Segmentation\\Patient 0\\'; % Escape backslashes
timeFrames = {'T_0','T_10','T_20', 'T_30', 'T_40','T_50','T_60','T_70','T_80', 'T_90'};

for t = 1:length(timeFrames)
    timeFrame = timeFrames{t};
    folderPath = fullfile(basePath, timeFrame, 'CT');
    filePattern = fullfile(folderPath, '*.dcm');
    dicomFiles = dir(filePattern);

    % ROI
    fullFileName = fullfile(folderPath, dicomFiles(1).name);
    dicomImage = dicomread(fullFileName);
    imshow(dicomImage, []);
    [croppedImage, rect] = imcrop;
    waitforbuttonpress;
    rect = round(rect);
    close all; 

    % Array 3D for memorizing the volume later
    vol = zeros(rect(4), rect(3), length(dicomFiles), 'logical'); 

    % Processing the DICOM files
    for k = 1:length(dicomFiles)
        fullFileName = fullfile(folderPath, dicomFiles(k).name);
        fprintf(1, 'Reading %s\n', fullFileName);
        dicomImage = dicomread(fullFileName);

  % Add Gaussian noise to the image
        noisy_image = imnoise(dicomImage,'gaussian',0,10^-4);

       
        croppedImage = noisy_image(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));
       
        % Edge Detection with Gradient mehod using Sobel Kernel 
        sobelX = fspecial('sobel');
        sobelY = sobelX';
        gradientImageX = imfilter(double(croppedImage), sobelX, 'replicate');
        gradientImageY = imfilter(double(croppedImage), sobelY, 'replicate');
        magnitude = sqrt(gradientImageX .^ 2 + gradientImageY .^ 2);
        edgeImage = magnitude > (0.1 * max(magnitude(:))); 

        
        edgeImage = imfill(edgeImage, 'holes');
        edgeImage = bwareaopen(edgeImage, 500);
        lungsMask = bwareaopen(edgeImage, 5000); 
        Lungs = logical(lungsMask); 

        % Plotting
        if k == 1
            figure;
            subplot(1, 4, 1);
            imshow(croppedImage);
            title('Cropped Image');

            subplot(1, 4, 3);
            imshow(~Lungs); 
            title('Cleaned Lungs');

            subplot(1, 4, 4);
            imhist(croppedImage);
            title('Histogram of Cropped Image');
        end

        % Applying the mask to each slice
        vol(:,:,k) = Lungs;
    end

    % 3D Volume
    if exist('volshow', 'file')
        figure;
        volshow(vol);
    end
end

