function [ regI , outT ] = non_rigid_core( I, R, inT, alpha )
%NON_RIGID_CORE Implements one step of the non rigid registration of 2 input matrices: 
% I is the image and 
% R is the "reference" image
% inT is the input Transformation matrix - set this to zeros(n) at the
% beginning
%alpha between 0 and 1 - sets aggressiveness
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
    return 0;
    msg = 'Trying to register unequal matrices'
    error(msg)%There is a big error!!
end



[Sx, Sy] = gradient(R); % Sx and Sy are the gradients of the Reference image
d2Sx = norm(Sx)*norm(Sx); %Square of determinant
d2Sy = norm(Sy)*norm(Sy);

thirion = 0; %based on Smart Align paper, thirion =1 would be more complete
%%%%Lines for Thirion method:
if thirion ==1
    d2Sx = d2Sx + norm(I-R)*norm(I-R);
    d2Sy = d2Sy + norm(I-R)*norm(I-R);
end
 
Sx = Sx/d2Sx; %normalize Sx - gradient of reference image R 
Sy = Sy/d2Sy; %normalize Sy - gradient of reference image R
%%%%Lines for Thirion method:

%%%%% End of Thirion method

  
[Mx, My] = gradient(I); % gradient of the "moving" image
d2Mx = norm(Mx)*norm(Mx);
d2My = norm(My)*norm(My);


%%%%Lines for Thirion method:
if thirion == 1
    d2Mx = d2Mx + norm(I-R)*norm(I-R);
    d2My = d2My + norm(I-R)*norm(I-R);
end
%%%%% End of Thirion method


Mx = Mx/d2Mx;  % normalized gradient of the moving image
My = My/d2My;
%%%%
%%%%End of future Thirion update section


dTx = - (I-R)*(Mx + Sx); %incremental transformation field for the moving image
dTy = - (I-R)*(My + Sy); 

%
%Section to smoothen the dTx and dTy matrices
%
smth_dTx = dTx; %smoothened transformation field.
smth_dTy = dTy;
%
%Need to implement later
%

outT = inT + alpha*()

end

