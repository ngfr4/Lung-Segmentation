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

        % The cropping based on the ROI previously chosen
        croppedImage = dicomImage(rect(2):(rect(2)+rect(4)-1), rect(1):(rect(1)+rect(3)-1));

        % Edge Detection with Gradient mehod using Sobel Kernel 
        sobelX = fspecial('sobel');
        sobelY = sobelX';
        gradientImageX = imfilter(double(croppedImage), sobelX, 'replicate');
        gradientImageY = imfilter(double(croppedImage), sobelY, 'replicate');
        magnitude = sqrt(gradientImageX .^ 2 + gradientImageY .^ 2);
        edgeImage = magnitude > (0.1 * max(magnitude(:)));

        % Filling the Holes
        edgeImage = imfill(edgeImage, 'holes');

        % Removing the artifacts, choosing the biggest objects
        edgeImage = bwareaopen(edgeImage, 500);
        lungsMask = bwareaopen(edgeImage, 5000); 

        % Mask to select ONLY the lungs
        Lungs = logical(lungsMask); 

        % Plotting
        if k == 1
            figure;
            sgtitle('Segmentation with Gradient+Sobel Kernel')
            subplot(1, 3, 1);
            imshow(croppedImage);
            title('Cropped Image');

            subplot(1, 3, 2);
            imshow(edgeImage); 
            title(' Lungs');

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
