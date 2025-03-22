%% 
% MANUAL THRESHOLDING WITHOUT THE HISTOGRAM EQUALIZATION
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

    % Cropping based on the ROI chosen from the first image
    croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));

    % Normalizing the image for the histogram then
    croppedImageDouble = double(croppedImage);
    normalizedCroppedImage = (croppedImageDouble - min(croppedImageDouble(:))) / (max(croppedImageDouble(:)) - min(croppedImageDouble(:)));

    % SEGMENTATION with Manual Thresholding (chosed based on the histogram)
    manualThreshold = 0.3; 
    Segmented_Lungs = imbinarize(normalizedCroppedImage, manualThreshold);
    
    %  we invert the binary image (Segmented_Lungs) to have lungs as white (1) and background as black (0).
    Segmented_Lungs = ~Segmented_Lungs;
    
     Lungs = bwareafilt(Segmented_Lungs, 2);
     Lungs = bwareaopen(Lungs, 500);
     Lungs = imfill(Lungs, 'holes');
     Lungs = bwareaopen(Lungs, 5000);

    % Applying the mask to each slice
    vol(:,:,k) = Lungs;
        % Plotting
    if k == 1

        figure
        imhist(normalizedCroppedImage);
        title('Histogram of Cropped Image');

        figure;
        subplot(1, 3, 1);
        imshow(croppedImage);
        title('Cropped Image');

        subplot(1, 3, 2);
        imshow(Segmented_Lungs);
        title('Segmented Lungs');

        subplot(1, 3, 3);
        imshow(Lungs); 
        title('Cleaned Lungs');

    end

      
    end

    % 3D Volume
    if exist('volshow', 'file')
        figure;
        volshow(vol);
    end
end
%%
% MANUAL THRESHOLDING WITH THE HISTOGRAM EQUALIZATION

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

        % Cropping based on the ROI chosen from the first image
        croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));

        % Histogram Equalization
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

        % Applying the mask to each slice
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
            imshow(Lungs); 
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

