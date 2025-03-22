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

        % The cropping based on the ROI previously chosen
        croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));
        
        % SEGMENTATION with Otsu's thresholding
        level = graythresh(croppedImage);
        Segmented_Lungs = imbinarize(croppedImage, level);

        %Selecting the main two black areas (the lungs without the borders)
        Lungs = bwareafilt(~Segmented_Lungs, 2);
        % Filling the holes in the lungs for the volume later
        Lungs = imfill(Lungs, 'holes');
        % Remove small black spots that are not part of the lungs
        Lungs = bwareaopen(Lungs, 500);  % Removes small objects (noise)
        % Keep only the largest connected components (the lungs)
        Lungs = bwareaopen(Lungs, 5000); % Keeps only large objects (lungs)

        % Plotting
    if k == 1
        figure;
        sgtitle('Segmentation with Otsu method')
        subplot(1, 3, 1);
        imshow(croppedImage);
        title('Cropped Image');

        subplot(1, 3, 2);
        imshow(Segmented_Lungs);
        title('Segmented Lungs');

        subplot(1, 3, 3);
        imshow(~Lungs); 
        title('Cleaned Lungs');

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

%% PLOTTING THE VOLUME IN TIME
clear all;
close all;
clc;

basePath = 'C:\Users\micol\OneDrive\Documenti\MATLAB\MI\PROJECT\Segmentation\Patient 0\';
timeFrames = {'T_0','T_10','T_20', 'T_30', 'T_40','T_50','T_60','T_70','T_80', 'T_90'};


totalLungVolumes = zeros(1, length(timeFrames));

firstFolderPath = fullfile(basePath, timeFrames{1}, 'CT');
firstFilePattern = fullfile(firstFolderPath, '*.dcm');
firstDicomFiles = dir(firstFilePattern);
firstFullFileName = fullfile(firstFolderPath, firstDicomFiles(1).name);
firstDicomImage = dicomread(firstFullFileName);
imshow(firstDicomImage, []);
[croppedImage, rect] = imcrop;
waitforbuttonpress;
rect = round(rect);
close all; 

for t = 1:length(timeFrames)
    timeFrame = timeFrames{t};
    folderPath = fullfile(basePath, timeFrame, 'CT');
    filePattern = fullfile(folderPath, '*.dcm');
    dicomFiles = dir(filePattern);

    % Array for storing cross-sectional areas
    sliceArea = zeros(1, length(dicomFiles)); 

    for k = 1:length(dicomFiles)
        fullFileName = fullfile(folderPath, dicomFiles(k).name);
        fprintf(1, 'Reading %s\n', fullFileName);
        dicomImage = dicomread(fullFileName);

        % Applying the same cropping based on the ROI previously chosen
        croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));
        
        % SEGMENTATION with Otsu's thresholding
        level = graythresh(croppedImage);
        Segmented_Lungs = imbinarize(croppedImage, level);
        Lungs = bwareafilt(~Segmented_Lungs, 2);
        Lungs = imfill(Lungs, 'holes');
        Lungs = bwareaopen(Lungs, 500);  
        Lungs = bwareaopen(Lungs, 5000); 

        % Calculate and store the cross-sectional area for this slice
        sliceArea(k) = bwarea(Lungs);
    end

    % Calculating total lung volume for this time frame
    totalLungVolumes(t) = sum(sliceArea);
end

% Plot the overall lung volume over time frames
figure;
plot(0:10:90, totalLungVolumes, '-o');
xlabel('Time Frame');
ylabel('Total Lung Volume');
title('Lung Volume Over Time');
%%


