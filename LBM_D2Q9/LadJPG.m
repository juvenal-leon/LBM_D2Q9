clear all
fullFileName = '/Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/mediosporosos/medio-1-crop-Scaled.jpg';
if exist(fullFileName, 'file')
  rgbImage = imread(fullFileName);
  
  binaryImage = im2bw(rgbImage, 0.3);
  subplot(1, 3, 1);
  imshow(rgbImage)
  
  subplot(1, 3, 2);
  imshow(binaryImage);
  
  binaryIM_Invert=imcomplement(binaryImage);
  subplot(1, 3, 3);
  imshow(binaryIM_Invert);
  
  dlmwrite('/Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/mediosporosos/matrix.txt',binaryIM_Invert,' ');
  
else
  errorMessage = sprintf('Image file does not exist:\n%s', fullFileName);
  uiwait(warndlg(errorMessage));
end