function [ regI  ] = non_rigid_core( I, R, alpha, a )
%NON_RIGID_CORE Implements one step of the non rigid registration of 2 input matrices: 
%
%I is the moving image and 
%
%R is the "reference" image
%
% alpha is a Thirion parameter
%
% a sets aggressivness, and should be between 0 and 1
%
%This function should go inside a loop until converged 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thirion paper described by Dirk-Jan Kroon and Slump allows you to make
% the matrix more stable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on "Smart Align" algorithm. 
% As described by Lewys Jones in : 10.1186/s40679-015-0008-4 along with
% descriptions in paper by Cachier on demons algorithm
% Written on March 14, 2016

%I = I/norm(I);
[nR, nC] = size(I); %nR and nC are number of rows and columns in the images
[xxx, yyy] = meshgrid(0:1.0/(nC-1):1, 0:1.0/(nR-1):1); 
x_I = xxx;
y_I = yyy;


if size(I) ~= size(R)
    %return 0;
    msg = 'Trying to register matrices with a size mismatch'
    error(msg)%There is a big error!!
end

[gRx, gRy] = gradient(R);
[gIx, gIy] = gradient(I);


imval_dif = I - R; %difference in intensity values of static and moving images at this pixel 
              
uX = imval_dif.*(gRx)./(norm(gRx)^2 + norm(gRy)^2 + (alpha^2)*(imval_dif.^2) )+ imval_dif.*gIx./(norm(gIx)^2 + norm(gIy)^2+ alpha^2*(imval_dif.^2) ); 
uY = imval_dif.*(gRy)./(norm(gRx)^2 + norm(gRy)^2 + (alpha^2)*(imval_dif.^2) )+ imval_dif.*gIy./(norm(gIx)^2 + norm(gIy)^2 + alpha^2*(imval_dif.^2) ); 
 
%Multiply shift by "aggressiveness" 
dX = a*uX;
dY = a*uY;

%Remove all NaN values from the dX, dY matrices
dX(isnan(dX)) = 0;
dY(isnan(dY)) = 0;

%Now maybe we want to smoothen dX and dY !!
smth_dX = dX;
smth_dY = dY;

%These are the (x,y) coordinates of each pixel of I after it has been moved
%to it's corresponding location in R
%each is of size [nR, nC]
x_I_reg = x_I + smth_dX;
y_I_reg = y_I + smth_dY;


%Moving image I matrix does not change because each pixel still has the same intensity
%now resample to produce I_reg

%data first converted to a N x 3 matrix, 
%where each row is one point, (xi,yi,zi)
%This is done so scatteredInterpolant class can be used

xx = reshape(x_I_reg, [nC*nR, 1]);
yy = reshape(y_I_reg, [nC*nR, 1]);
zz = reshape(I, [nC*nR, 1]);

%xx, yy and zz are row vectors containing the (X, Y and Z) positions of
%each displaced pixel. 
%adjust this variable to choose the method of interpolation! You can either
%set it to 'griddata' or to 'scattin'
resamp = 'scattin';

zzz = 0*xxx; % This array will store new image values- same dimension as I

%Pick resampling strategy here
if strcmp(resamp, 'griddata') % using grid data matlab method
    zzz = griddata(xx,yy,zz,xxx, yyy);
elseif strcmp(resamp, 'scattin') % using scattered interpolant matlab method - this is recommended
    F = scatteredInterpolant(xx,yy,zz, 'natural'); % this is the interpolant function
    zzz = F(xxx,yyy);
end

%Now convert (X,Y,Z) data to image, regI
%Is this step needed??
regI = reshape(zzz , [nR,nC]);

% 
% figure
% subplot(1,3,1)
% imagesc(R); title('Orignial Reference Image'); colormap gray;
% subplot(1,3,2);
% imagesc(I); title('Original Moving Image');colormap gray;
% subplot(1,3,3);
% imagesc(regI); colormap gray;

end

