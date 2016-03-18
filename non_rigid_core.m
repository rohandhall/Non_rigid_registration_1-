function [ regI  ] = non_rigid_core( I, R, alpha, a )
%NON_RIGID_CORE Implements one step of the non rigid registration of 2 input matrices: 
%
%I is the moving image and 
%
%R is the "reference" image
%
% inT is the input Transformation matrix - set this to zeros(n) at the
% beginning
% alpha is a Thirion parameter
%
% a sets aggressivness, and should be between 0 and 1
%This function should go inside a loop until converged 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thirion paper described by Dirk-Jan Kroon and Slump allows you to make
% the matrix more stable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on "Smart Align" algorithm. 
% As described by Lewys Jones in : 10.1186/s40679-015-0008-4 along with
% descriptions in paper by Cachier on demons algorithm
% Written on March 14, 2016

[nR, nC] = size(I); %nR and nC are number of rows and columns in the images
if size(I) ~= size(R)
    %return 0;
    msg = 'Trying to register matrices with a size mismatch'
    error(msg)%There is a big error!!
end
%Create matrices x_I and y_I which store the x and y coordinates of each
%pixel. These matrices are the same size as image I
for r= 1:nR
    for c= 1:nC
        y_I(r,c) = (r-1)/(nR-1);
        x_I(r,c) = (c-1)/(nC-1);
    end
end

% [Sx, Sy] = gradient(R); % Sx and Sy are the gradients of the Reference image. Each are the same size as the image I
% [Mx, My] = gradient(I); % gradient of the "moving" image

%I didnt like the names S and M... changed them to gRx and gIx (Gradient of
%R or Gradient of I )
[gRx, gRy] = gradient(R);
[gIx, gIy] = gradient(I);

version = 2; %Realized there was a bug in version 1.

thirion = 1; 


if version ==2
    
    if thirion == 1
        for r = 1: nR
            for c = 1:nC
                gradR = [gRx(r,c), gRy(r,c)]; %gradient of reference image at this pixel 
                gradI = [gIx(r,c), gIy(r,c)]; %gradient of moving image at this pixel
                imval_dif = I(r,c) - R(r,c); %difference in intensity values of static and moving images at this pixel   
                u{r,c} = imval_dif*gradR/(norm(gradR)^2 + (alpha^2)*(imval_dif^2) )+ imval_dif*gradI/(norm(gradI)^2 + alpha^2*(imval_dif^2) ); 
            end
        end
    else
        for r = 1: nR
            for c = 1:nC
                gradR = [gRx(r,c), gRy(r,c)]; %gradient of reference image at this pixel
                gradI = [gIx(r,c), gIy(r,c)]; %gradient of moving image at this pixel
                imval_dif = I(r,c) - R(r,c); %difference in intensity values of static and moving images at this pixel   
                u{r,c} = imval_dif*gradR/(norm(gradR)^2)+ imval_dif*gradI/(norm(gradI)^2); 
            end
        end
    
    end
    
end


%Version 2: march 18, 2016 :
%u{r,c} contains the displacement of pixel {r, c}
%each cell is a vector [ dX, dY ]  
dX = 0*I;
dY = 0*I;

for r=1:nR
    for c = 1:nC
        tempdX = a*u{r,c}(1);
        tempdY = a*u{r,c}(2);
        if isnan(tempdX)
            tempdX  = 0;
        end
        if isnan(tempdY)
            tempdY =0;
        end
        
        dX(r,c) = tempdX;
        dY(r,c) = tempdY;
    end
end

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
rr =1;
for r = 1:nR
    for c = 1:nC
        xx(rr) = x_I_reg(r,c);
        yy(rr) = y_I_reg(r,c);
        zz(rr) = I(r,c);
    end
end

%NEED TO WORK THIS OUT !!
F = scatteredInterpolant(xx,yy,zz, 'natural'); % this is the interpolant function
regI = 0*I;
for r=1:nR
    for c=1:nC
        regI(r,c) = F(x_I_reg(r,c) , y_I_reg(r,c));
    end
end

end

