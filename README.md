# FSC_SMLM package
Fourier shell correlation from XYZ coordinates

# Requirement
Matlab 2019 or later

Dipimage toolbox 2.9 ([https://diplib.org/](https://diplib.org/download_2.9.html))

FIREfunctions (https://www.nature.com/articles/nmeth.2448)


# Demo_FSC.m
Demo code for Fourier sell correlation (FSC) analysis
CSV localization file is selected throuhg UI file selection. In CSV file,
the 3-5 columns in CSV file must be X, Y, and Z coordinates in nm. FSC is
perfomed though several subregions to accelerate analysis and save
memory. To prevent biased estimation, FSC is performed only if
localizatyion number inside targeted region is larger than 200.

# FSC_SMLM.m
Input :

      X,Y,Z: Nx1 position array of coordinates in specific axis. Unit in nm. 
      Np: Isotropic size of subvolume for FSC. Unit in Pixel.
      Para: struct of parameters related to camera pixel size.
          sz: FOV size in Pixel
          psz: Effective pixel size in nm/pixel
          zm: Zoom-in ratio for FSC.

Output :

      S: struct of output results
          fsc: Mx1 double array of FSC correlation result. Contain only frequencies lower than Nyquist frequency.
          fsc_s: Mx1 double array of fitted FSC correlation result with fitting type of smoothingspline with sampling frequency of 1/3 data. Contain only frequencies lower than Nyquist frequency.
          xu: Mx1 double array of corresponding FSC frequency in 1/pixel. Contain only frequencies lower.
          nu: Mx1 double array of corresponding FSC frequency in 1/nm. Contain only frequencies lower.
          r_17: 1/7-threshold resolution of FSC result in nm.

Fourier Shell Correlation (FSC) is calculatedby two given XYZ coordinates. 
Dataset will be equally split into two subdataset, rendered as 3D images,
and applied 3D Tukey window before performing FSC.

Created by HAO-CHENG GAO

Purdue University, 610 Purdue Mall, West Lafayette, IN 47907

gao561@purdue.edu

Copyright. June 2, 2023.


# References:
Nieuwenhuizen, R., Lidke, K., Bates, M. et al. "Measuring image resolution in optical nanoscopy" Nat Methods 10, 557–562 (2013).

Huang, F., Sirinakis, G., et al. "Ultra-High Resolution 3D Imaging of Whole Cells" Cell 166, 4, 1028-1040 (2016)
      
      
