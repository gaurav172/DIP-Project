% Code for Illuminant Spectra-based Source Separation Using Flash Photography
% This code is based on the algorithm proposed in the paper

% "Illuminant Spectra-based Source Separation Using Flash Photographye", CVPR 2018
% Zhuo Hui, Kalyan Sunkavalli, Sunil Hadap, Aswin C. Sankaranarayanan


% When you use the code to build your algorithm, please cite this paper. 
% 
% Please contact the author Zhuo Hui if you have any problems with the code
% huizhuo1987@gmail.com
% 
% Copy rights reserved by the authors listed above.

%% This function is to implement the esimtation of lighting coefficients

% Input
% 1. L: measured illuminant spectra data N by M, 
%    N: the samples for the spectra, M the number of samples in the database   
% 2. R: measured reflectance spectra data N by M
% 3. C: camera response spectra
% 4. option: the difference methods in generating the basis

%   

% Output
% 1. UR, UL: the generated basis for reflectance and illumination,
%    respectively



function [UR, UL] = genBase(L, C, R, option)

switch option
    case 'joint'
        t1 = diag(C*ones(3,1))
        t2 = R' * t1 * L
        [VR, VK, VL] = svds(t2, 3);
        UR = orth(R*VR);
        UL = orth(L*VL);
          
    case 'pca'

        [UR, SigmaR, VR] = svds(R, 3);
        
        [UL, SigmaL, VL] = svds(L, 3);
        UR = orth(UR);
        UL = orth(UL);
    case 'wpca'
        a1 = R.*repmat((C(:, 1)),1,size(R,2));
        a2 = R.*repmat((C(:, 2)),1,size(R,2));
        a3 = R.*repmat((C(:, 3)),1,size(R,2));
        a4 = L.*repmat((C(:, 1)),1,size(L,2));
        a5 = L.*repmat((C(:, 2)),1,size(L,2));
        a6 = L.*repmat((C(:, 3)),1,size(L,2));
        R_total = [a1,a2,a3];
        L_total = [a4,a5,a6];
        [UR, UK, UP] = svds(R_total, 3);       
        [UL, UK, UP] = svds(L_total, 3);   

        
        UR = orth(UR);
        UL = orth(UL);
    case 'wpca1'
        a1 = R.*repmat((C(:, 1)), 1, size(R, 2));
        a2 = R.*repmat((C(:, 2)), 1, size(R, 2));
        a3 = R.*repmat((C(:, 3)), 1, size(R, 2));
        a4 = L.*repmat((C(:, 1)), 1, size(L, 2));
        a5 = L.*repmat((C(:, 2)), 1, size(L, 2));
        a6 = L.*repmat((C(:, 3)), 1, size(L, 2));

        [UR(:, 1),U11,U12] = svds(a1,1);
        [UR(:, 2),U21,U22] = svds(a2,1);
        [UR(:, 3),U31,U32] = svds(a3,1);

        [UL(:, 1),U41,U42] = svds(a4,1);
        [UL(:, 2),U51,U52] = svds(a5,1);
        [UL(:, 3),U61,U62] = svds(a6,1);
        
        UR = orth(UR);
        UL = orth(UL);
end


end