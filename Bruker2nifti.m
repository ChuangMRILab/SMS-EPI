function Bruker2nifti(pathData,pathDest,crop,scale)
% Convert Bruker data to nifti file
%
% Input:
%       pathData -  path to scan folder
%       pathDest -  path to put the nifti file, default is current path
%       crop     -  set to 1 to perform initial cropping of the image
%       scale    -  voxel size scaling factor, default is 10
%
     

if nargin < 4
    scale = 10;
end

if nargin < 3
    crop = 0;
end

if nargin < 2
    pathDest = pwd;
end

% Build nifti struct from scan folder
[data,nifti] = Bruker2nifti_smsEPI(pathData,scale,crop);

% Parse file name
[pathstr,Enum,fileName] = fileparts(pathData);
newFileName             = [pathDest filesep data.info.descrip 'X' Enum 'P1.nii.gz'];

% Save nifti to file
save_untouch_nii(nifti,newFileName);
