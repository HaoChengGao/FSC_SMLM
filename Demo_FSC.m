% Demo code for Fourier sell correlation (FSC) analysis
% CSV localization file is selected throuhg UI file selection. In CSV file,
% the 3-5 columns in CSV file must be X, Y, and Z coordinates in nm. FSC is
% perfomed though several subregions to accelerate analysis and save
% memory. To prevent biased estimation, FSC is performed only if
% localizatyion number inside targeted region is larger than 200.

% Last modified on June 2, 2023 by HAO-CHENG GAO
%=========================================================================

addpath('\Function')
addpath('\Function\FIREfunctions')

[file, path] = uigetfile('C:\*.csv','Select localization csv file:','MultiSelect','off');


a = readmatrix(fullfile(path,file));
obj = struct('x',a(:,2),'y',a(:,3),'z',a(:,4));

%%  Parameter

para.sz = 297;      % FOV size in Pixel
para.psz = 159;     % Effective pixel size in nm/pixel
para.zm = 10;       % Zoom-in ratio for FSC. The value will influence Nyquist frequency(Nf),e.g. Nf = para.zm/para.psz/2 . Default as 10.

N = 11;             % Subregion number for analysis. Boundary subregions will not being used in analysis. Default as 11;


%%
Lm = linspace(-para.sz*para.psz/2, para.sz*para.psz/2, N+1);
Np = ceil(diff(Lm(1:2))/para.psz)*para.zm;

FSC_list = cell(N-2,N-2);


for i = 6%2:N-1
    for j = 6%2:N-1
        % Select localizations inside targeted subregion.
        msk = obj.x>=Lm(i) & obj.x<Lm(i+1) & obj.y>=Lm(j) & obj.y<Lm(j+1);
        Nij = sum(msk(:));
        disp([ 'Subregion (' num2str(i) ',' num2str(j) ') with ' num2str(Nij) ' localizations']);

        % Perform FSC only if localization number inside targeted subregion
        % is larger than 200 (preventing biased estimation)
        if (Nij>200)
            fsc = FSC_SMLM(obj.x(msk), obj.y(msk), obj.z(msk), Np, para);
            FSC_list(j-1,i-1) = {fsc};
        end
    end
end



