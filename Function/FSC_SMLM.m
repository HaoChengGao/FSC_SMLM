function S = FSC_SMLM(X, Y, Z, Np, para)
% Input :
%
%       X,Y,Z: Nx1 position array of coordinates in specific axis. Unit in
%       nm. 
%       Np: Isotropic size of subvolume for FSC. Unit in Pixel.
%       Para: struct of parameters related to camera pixel size.
%           sz: FOV size in Pixel
%           psz: Effective pixel size in nm/pixel
%           zm: Zoom-in ratio for FSC.
%
% Output :
%       S: struct of output results
%           fsc: Mx1 double array of FSC correlation result. Contain only 
%                frequencies lower than Nyquist frequency.
%           fsc_s: Mx1 double array of fitted FSC correlation result with 
%                  fitting type of smoothingspline with sampling frequency 
%                  of 1/3 data. Contain only frequencies lower than Nyquist 
%                  frequency.
%           xu: Mx1 double array of corresponding FSC frequency in 1/pixel.
%               Contain only frequencies lower.
%           nu: Mx1 double array of corresponding FSC frequency in 1/nm.
%               Contain only frequencies lower.
%           r_17: 1/7-threshold resolution of FSC result in nm.
%
% Fourier Shell Correlation (FSC) is calculatedby two given 3d-coordinates. 
% Dataset will be equally split into two subdataset, rendered as 3D images,
% and applied 3D Tukey window before performing FSC.

% Created by HAO-CHENG GAO
%           Purdue University 
%           610 Purdue Mall
%           West Lafayette
%           IN 47907
%           gao561@purdue.edu
%
% Copyright. June 2, 2023.
%
% References:
% Nieuwenhuizen, R., Lidke, K., Bates, M. et al. "Measuring image
%       resolution in optical nanoscopy" Nat Methods 10, 557–562 (2013).
% F. Huang, G. Sirinakis, et al. "Ultra-High Resolution 3D Imaging of Whole
%       Cells" Cell 166, 4, 1028-1040 (2016)
%=========================================================================

% Split dataset into two and align them to coordinate origin
indx1 = mod(1:length(X),2)==1;
indx2 = ~indx1;

x1 = X(indx1)-min(X)+1;
x2 = X(indx2)-min(X)+1;

y1 = Y(indx1)-min(Y)+1;
y2 = Y(indx2)-min(Y)+1;

z1 = Z(indx1)-min(Z)+1;
z2 = Z(indx2)-min(Z)+1;

xm1 = x1/para.psz;
xm2 = x2/para.psz;
ym1 = y1/para.psz;
ym2 = y2/para.psz;
zm1 = z1/para.psz;
zm2 = z2/para.psz;

%%



% Turkey Window
xt = 1:Np;
[Xt,Yt,Zt] = meshgrid(xt,xt,xt);
Xt1 = sin(4*pi*Xt/Np).^2;
Xt1(:,ceil(Np/8):floor(Np*7/8),:)=1;
Yt1 = sin(4*pi*Yt/Np).^2;
Yt1(ceil(Np/8):floor(Np*7/8),:,:)=1;
Zt1 = sin(4*pi*Zt/Np).^2;
Zt1(:,:,ceil(Np/8):floor(Np*7/8))=1;
WindowTK = Xt1.*Yt1.*Zt1;
        
imstr1 = binlocalizations3D([xm1,ym1,zm1],Np,Np,Np,para.zm).*WindowTK;
imstr2 = binlocalizations3D([xm2,ym2,zm2],Np,Np,Np,para.zm).*WindowTK;

%%  FSC

% Perform FSC and fit the result with smoothing spline curve
S = FSC(double(imstr1),double(imstr2));

nu = S.xu/para.psz*para.zm;
sp =round(linspace(1,length(nu),round(length(nu)/3))); 
ff =fit(nu(sp), S.fsc(sp),'smoothingspline'); 
fsc_s = feval(ff,nu);
Rn2 = find(fsc_s<1/7,1,'first');
S.r_17 = 1/nu(Rn2);
S.fsc_s = fsc_s;
S.nu = nu;

% Figure to plot FSC result with fitted FSC curve and 1/7-threshold resolution
figure
hold on
plot(nu, S.fsc,'.','Markersize',15)
plot(nu, fsc_s,'linew',3)
line([nu(Rn2) nu(Rn2)],[0 1],'linestyle','-','color',[0.6 0.6 0.6],'linew',3)
text(nu(Rn2)+0.001,1/6,[num2str(1/nu(Rn2),'%2.1f') ' nm'],'HorizontalAlignment','left','VerticalAlignment','baseline')
xlim([0 max(nu)])
xlabel('Spatial frequency (1/nm)')
ylabel('Correlation')
legend('FSC','FSC_s_m_o_o_t_h_s_p_l_i_n_e','1/7 threshold')


end

function S = FSC( im1, im2)

% Fourier transform both satasets
f1 = fftn(im1);
f2 = fftn(im2);
N = size(f1,1);

% generate shell label to Fourier domain 
xi = 0;
for jj = 1:3
    xi = bsxfun( @plus, xi, reshape((-floor(0.5*N):ceil(0.5*N)-1).^2, [ones([1,jj-1]), N, ones([1,2-jj])]));
end
s = round(sqrt(ifftshift(xi)))+1;
Ns = [max(s(:)), 1];

% Compute correlation on shells considering limitation of Nyquist frequency
S.fsc =    real(accumarray(s(:), f1(:).*conj(f2(:)), Ns)) ./ sqrt(accumarray(s(:), abs(f1(:)).^2, Ns) .* accumarray(s(:), abs(f2(:)).^2, Ns));
S.fsc = S.fsc(1:floor(ceil((N+1)/2)));
S.xu = (0:numel(S.fsc)-1)'/numel(S.fsc)/2;


end
