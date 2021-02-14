# Flake Dimensional Analysis GUI
Flake Identification and measurement GUI using Matlab's image processing toolbox to find and quantify flakes in an AFM .IBW scan.
### The Image Processing Toolbox is required

## Quick Install (MATLAB R2012b? and Beyond)
Under ".\Flake Analysis App" open "AFM Image Analysis" using MATLAB, this will install the GUI. 
The app can later be accessed within MATLAB by going to the APPS tab and selecting "AFM Image Analysis".

## Manual Install and edits
To read .IBW files, this GUI relies on Jakub Bialek's [file converter](https://www.mathworks.com/matlabcentral/fileexchange/42679-igor-pro-file-format-ibw-to-matlab-variable). 
The APP already includes these as part of the package, but to successfully run the .mlapp you will need to download these files.
1. Download all IBW parsing files from the link above
2. Put all .m files within the directory with a folder called "IBW"


Note that older versions of MATLAB will have trouble opening .mlapp programs, used in this project to design to GUI. This [thread](https://www.mathworks.com/matlabcentral/answers/404815-how-to-get-code-from-a-mlapp-file-using-an-earlier-matlab-version)
contains some work arounds. The image processing and logic components of this project can be called exclusively through the "getAFMFlakes" function. The arguments for this function are as follows...
```
%getAFMFlakes 
% Returns a structure with various flake data from an AFM file
%
% ARGUMENTS
%   file            File path to IBW AFM scan
%   contrast        Toggle for increasing contrast during image processing
%   edgeMethod      Algorithm to detect edges in image for flake identification
%   fudgeFactor     Sensitivity of edge detection
%   filterRadius    Noise removal filter radius (Weiner filter)
%   dilateRadius    Dilation factor of identified flakes
%
% RETURNS
% afmImgInfo.   area            List of flake areas
%               length          List of flake lengths
%               meanThick       List of mean flake thicknesses
%               medianThick     List of median flake thicknesses
%               maxThick        List of maximum flake thicknesses
%               volume          List of flake volumes
%               rawAFMData      The raw height information from AFM
%               labeledFlakes   An image with identified flakes outlined 
%               labeledFlakesNum An image with identified flakes colored and numbered
%               flakeInfo       A 2D array of sorted flake information for GUI display
```
