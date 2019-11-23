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


%% This function is to implement the lighting separation

% Input
% 1. im_nf: no flash image M*N*3
% 2. im_f: flash image M*N*3
% 3. mask: mask for the image M*N
% 4. opt includes: 
%    1. UR: the basis of reflectance
%    2. UL: the basis of illumination
%    3. E: the matrix E as descirbed in our paper
%    4. f: the lighting coefficients for the flash lighting 
%    optional 
%    5. pMask: the mask for the pixels with purely colors
%    6. lambda: the weights for the ridge regression solver for beta
%    7. shadow_mask: optional. to mask the flash shdows
%    8. cutoff: the threshold values in thresholding the histogram

%   

% Output
% The output is in the structure format, which includes
% 1. alpha: the alpha image as described in our paper
% 2. beta: the beta image as described in our paper
% 3. im1, im2: separated images
% 4. im1_wb, im2_wb: white balanced images
% 5. illum1, illum2: estimated lighing coefficients
% 6. coeff: estimated relative shading

function results = solve_light_sep(im_nf, im_f, mask, opt)

%% initialize the variables
pMask = mask;
% if ~isfield(opt, 'pMask')
%     pMask = mask;
% else
if isfield(opt,'pMask')
    pMask = opt.pMask;
end

% if ~isfield(opt, 'lambda')
lambda = 1e-3;
% else
if isfield(opt,'lambda')
    lambda = opt.lambda;
end

cfactor = 1;

% if ~isfield(opt, 'cutoff')
%     cfactor = 1;
% else
if isfield(opt,'cutoff')
    cfactor = opt.cutoff;
end

% if ~isfield(opt, 'E')
%     error('Not found E matrix')
%     return
% else
E = opt.E;
% end

% if ~isfield(opt, 'f')
%     error('Not found flash coefficients')
%     return
% else
f = opt.f;
% end

% if ~isfield(opt, 'light_number')
%     error('Please specificy the number of lights in the scene')
%     return
% else
Q = opt.light_number;
% end

siz = size(im_nf);

%% generate difference image
D_img = im_f - im_nf;
D_img(D_img<0) = 0;
Temp_Mat = repmat(mask,[1 1 3]);
diff_img = D_img.*Temp_Mat;
% Diff_img = max(im_f - im_nf, 0).*repmat(mask,[1 1 3]);
% diff_img = diff_img/1000;
for i=1:3
    Amat(:,i) = E(:,:,i)*f;
end
% Amat(:, 1) = E(:,:,1)*f;
% Amat(:, 2) = E(:,:,2)*f;
% Amat(:, 3) = E(:,:,3)*f;

diff_img_vec = reshape(diff_img, [], 3);
diff_img_vec = diff_img_vec';
img_nf_vec = reshape(im_nf, [], 3)';
img_nf_vec = img_nf_vec;

%% get alpha image
alpha = (Amat');
alpha = alpha\diff_img_vec;

if ~isfield(opt, 'shadow_mask')
    alpha = shadowRM(im_nf, im_f, mask, reshape(alpha', siz), opt.shadow_mask)';
%     alpha = alpha';
end

trp = alpha';
results.alpha = reshape(trp, siz);
tx = sum(alpha.^2,1);
tx = sqrt(tx);
norm_alpha = tx;
tx = tx';
results.norm_alpha = reshape(tx, siz(1:2));

%% 
for kk = 1:size(diff_img_vec, 2)
    if mod(kk, 10000) == 1
        fprintf('\b\b\b\b\b\b\b%1.3f \n', 3*kk/prod(size(diff_img_vec)));
    end
    if (mask(kk))
        %% solve for the beta image
        for i=1:3
            Bmat(i,:) = alpha(:,kk)'*E(:,:, i);
        end
%         Bmat(1, :) = alpha(:, kk)'*E(:,:, 1);
%         Bmat(2, :) = alpha(:, kk)'*E(:,:, 2);
%         Bmat(3, :) = alpha(:, kk)'*E(:,:, 3);
        Bmat = Bmat/(1e-10+norm_alpha(kk));
        
        
        nf_intensity = img_nf_vec(:, kk);
        num_1 = Bmat'*Bmat;
        num_2 = lambda*eye(3);
        deno = Bmat'*nf_intensity;
        z = 0;
        if ~pMask(kk)
            z = 5;
%             beta(:, kk) = (Bmat'*Bmat+5*lambda*eye(3))\(Bmat'*nf_intensity);
        else
            z = 1;
%             beta(:, kk) = (Bmat'*Bmat+5*lambda*eye(3))\(Bmat'*nf_intensity);
        end
        beta(:,kk) = (num_1+z*num_2)\deno;
        if sum(isnan(beta(:, kk)))
            beta(:, kk) = zeros(3, 1);
        end
        beta_norm(kk) = norm(beta(:,kk));  
        gamma(:, kk) = beta(:,kk)/(1e-10+beta_norm(kk));
    else
        beta(:, kk) = zeros(3, 1);
        beta_norm(kk) = 0;
        gamma(:, kk) = zeros(3, 1);
    end
end


%% filter out the noise for gamma image
im_nf_t = im_nf;
im_nf_t = im_nf_t.*repmat(mask, [1 1 3]);
s_img = reshape(beta', siz);
smoothness = 0.0000001; %0.0001;%
% for i = 1:3
%     s_img(:,:,i) = imguidedfilter(s_img(:,:,i), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
% end
s_img(:,:,1) = imguidedfilter(s_img(:,:,1), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img(:,:,2) = imguidedfilter(s_img(:,:,2), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img(:,:,3) = imguidedfilter(s_img(:,:,3), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img = reshape(s_img, [], 3)';
% s_img = s_img';
sums = sum(s_img.^2,1);
sqsum = sqrt(sums);
beta_norm = sqsum;
deno = repmat(beta_norm,[3,1]);
deno = deno + 1e-10;
gamma = s_img./deno;

%%
results.beta = reshape(beta', siz);
results.beta_norm = reshape(beta_norm', siz(1:2));
results.gamma = reshape(gamma', siz);


%% estimated the illumination coefficients
mask_t = mask & pMask;
if Q == 2
    %% Two lights case using ransac
    [illum1, illum2] = est_two_light_coeff(gamma, mask_t, cfactor);
    illum = [illum1,illum2]
    L_mat = illum;
    %% reproject the shadings
    mR1 = (alpha'*E(:,:,1)*illum(:,1))'.*(beta_norm./(1e-10+norm_alpha));
    mR2 = (alpha'*E(:,:,1)*illum(:,2))'.*(beta_norm./(1e-10+norm_alpha));

    mG1 = (alpha'*E(:,:,2)*illum(:,1))'.*(beta_norm./(1e-10+norm_alpha));
    mG2 = (alpha'*E(:,:,2)*illum(:,2))'.*(beta_norm./(1e-10+norm_alpha));

    mB1 = (alpha'*E(:,:,3)*illum(:,1))'.*(beta_norm./(1e-10+norm_alpha));
    mB2 = (alpha'*E(:,:,3)*illum(:,2))'.*(beta_norm./(1e-10+norm_alpha));

    for kk = 1:size(img_nf_vec,2)
        Amat_ = [mR1(kk) mR2(kk); mG1(kk) mG2(kk); mB1(kk) mB2(kk)];
        bvec_ = img_nf_vec(:, kk);
        coeff(:, kk) = Amat_\bvec_;
    end

%%
    for i=1:2
        corr1 = coeff(1, :).*beta_norm./(1e-10+norm_alpha);
        %(coeff1./gamma_pred_norm).*beta_norm./(1e-10+norm_alpha);
        corr2 = coeff(2, :).*beta_norm./(1e-10+norm_alpha);
%         size(corr1)
        corr(i,:) = coeff(i,:).*beta_norm./(1e-10+norm_alpha);
%         size(corr(i,:))
    end
%     size(corr1)
%     size(corr(:,1))
    %(coeff2./gamma_pred_norm).*beta_norm./(1e-10+norm_alpha);
%     for i=1:2
%           MR(i,:)=(alpha'*E(:,:,1))'.*repmat(corr(i,:),[3 1]);
%           MG(i,:)=(alpha'*E(:,:,2))'.*repmat(corr(i,:),[3 1]);
%           MB(i,:)=(alpha'*E(:,:,3))'.*repmat(corr(i,:),[3 1]);
%     end
    MR1 = (alpha'*E(:,:,1))'.*repmat(corr(1,:), [3 1]);
    MR2 = (alpha'*E(:,:,1))'.*repmat(corr(2,:), [3 1]);
%     MR = [MR1,MR2];
%     size(MR)
    MG1 = (alpha'*E(:,:,2))'.*repmat(corr(1,:), [3 1]);
    MG2 = (alpha'*E(:,:,2))'.*repmat(corr(2,:), [3 1]);

    MB1 = (alpha'*E(:,:,3))'.*repmat(corr(1,:), [3 1]);
    MB2 = (alpha'*E(:,:,3))'.*repmat(corr(2,:), [3 1]);

    b= [MR1' MR2'; MG1' MG2'; MB1' MB2']\(reshape(img_nf_vec', [], 1));
    illum_1 = b(1:3);
    illum_2 = b(4:6);
    Illum_n = [b(1:3),b(4:6)];
    illum_1 = illum_1/norm(illum_1);
    illum_2 = illum_2/norm(illum_2);
    Illum_n(:,1) = Illum_n(:,1)/norm(Illum_n(:,1));
    Illum_n(:,2) = Illum_n(:,2)/norm(Illum_n(:,2));
    im1_new = [];
    im1_new(1, :) = (alpha'*E(:,:,1)*Illum_n(:,1))'.*corr(1,:);
    im1_new(2, :) = (alpha'*E(:,:,2)*Illum_n(:,1))'.*corr(1,:);
    im1_new(3, :) = (alpha'*E(:,:,3)*Illum_n(:,1))'.*corr(1,:);
%     for kk=1:3
%         im1_new(kk, :) = (alpha'*E(:,:,kk)*Illum_n(:,1))'.*corr1;
%     end
    im1_new = reshape(im1_new', size(im_nf));

    im2_new = [];
    im2_new(1, :) = (alpha'*E(:,:,1)*Illum_n(:,2))'.*corr(2,:);
    im2_new(2, :) = (alpha'*E(:,:,2)*Illum_n(:,2))'.*corr(2,:);
    im2_new(3, :) = (alpha'*E(:,:,3)*Illum_n(:,2))'.*corr(2,:);
%     for kk=1:3
%         im2_new(kk, :) = (alpha'*E(:,:,kk)*Illum_n(:,2))'.*corr2;
%     end
    im2_new = reshape(im2_new', size(im_nf));


    im1_wb = [];
    im1_wb(1, :) = (alpha'*E(:,:,1)*f)'.*corr(1,:);
    im1_wb(2, :) = (alpha'*E(:,:,2)*f)'.*corr(1,:);
    im1_wb(3, :) = (alpha'*E(:,:,3)*f)'.*corr(1,:);
%     for kk=1:3
%         im1_wb(kk, :) = (alpha'*E(:,:,kk)*f)'.*corr1;
%     end
    im1_wb = reshape(im1_wb', size(im_nf));

    im2_wb = [];
    im2_wb(1, :) = (alpha'*E(:,:,1)*f)'.*corr(2,:);
    im2_wb(2, :) = (alpha'*E(:,:,2)*f)'.*corr(2,:);
    im2_wb(3, :) = (alpha'*E(:,:,3)*f)'.*corr(2,:);
%     for kk=1:3
%         im2_wb(kk, :) = (alpha'*E(:,:,kk)*f)'.*corr2;
%     end
    im2_wb = reshape(im2_wb', size(im_nf));  

    results.im1 = im1_new;
    results.im2 = im2_new;   
    results.im1_wb = im1_wb;
    results.im2_wb = im2_wb;
    results.coeff = coeff;
%% Three lights case using ransac
else 
    [illum1, illum2, illum3] = est_three_light_coeff(gamma, mask_t, cfactor);
    L_mat = [illum1 illum2 illum3];
    coeff = .5*ones(Q, size(gamma, 2));

    coeff(:, mask > 0) = L_mat\gamma(:, mask > 0);
    
    illum1 = L_mat(:, 1); illum2 = L_mat(:, 2); illum3 = L_mat(:, 3);
    mR1 = (alpha'*E(:,:,1)*illum1)'.*(beta_norm./(1e-10+norm_alpha));
    mR2 = (alpha'*E(:,:,1)*illum2)'.*(beta_norm./(1e-10+norm_alpha));
    mR3 = (alpha'*E(:,:,1)*illum3)'.*(beta_norm./(1e-10+norm_alpha));
    
    mG1 = (alpha'*E(:,:,2)*illum1)'.*(beta_norm./(1e-10+norm_alpha));
    mG2 = (alpha'*E(:,:,2)*illum2)'.*(beta_norm./(1e-10+norm_alpha));
    mG3 = (alpha'*E(:,:,2)*illum3)'.*(beta_norm./(1e-10+norm_alpha));
    
    mB1 = (alpha'*E(:,:,3)*illum1)'.*(beta_norm./(1e-10+norm_alpha));
    mB2 = (alpha'*E(:,:,3)*illum2)'.*(beta_norm./(1e-10+norm_alpha));
    mB3 = (alpha'*E(:,:,3)*illum3)'.*(beta_norm./(1e-10+norm_alpha));
    
    for kk = 1:size(img_nf_vec,2)
        Amat_ = [mR1(kk) mR2(kk) mR3(kk); mG1(kk) mG2(kk) mG3(kk); mB1(kk) mB2(kk) mB3(kk)];
        bvec_ = img_nf_vec(:, kk);
        coeff(:, kk) = Amat_\bvec_;
    end    
        
    %%
    corr1 = coeff(1, :).*beta_norm./(1e-10+norm_alpha);
    corr2 = coeff(2, :).*beta_norm./(1e-10+norm_alpha);
    corr3 = coeff(3, :).*beta_norm./(1e-10+norm_alpha);
    MR1 = (alpha'*E(:,:,1))'.*repmat(corr1, [3 1]);
    MR2 = (alpha'*E(:,:,1))'.*repmat(corr2, [3 1]);
    MR3 = (alpha'*E(:,:,1))'.*repmat(corr3, [3 1]);
    
    MG1 = (alpha'*E(:,:,2))'.*repmat(corr1, [3 1]);
    MG2 = (alpha'*E(:,:,2))'.*repmat(corr2, [3 1]);
    MG3 = (alpha'*E(:,:,2))'.*repmat(corr3, [3 1]);
    
    MB1 = (alpha'*E(:,:,3))'.*repmat(corr1, [3 1]);
    MB2 = (alpha'*E(:,:,3))'.*repmat(corr2, [3 1]);
    MB3 = (alpha'*E(:,:,3))'.*repmat(corr3, [3 1]);
    
    b= [MR1' MR2' MR3'; MG1' MG2' MG3'; MB1' MB2' MB3']\(reshape(img_nf_vec', [], 1));
    illum_1 = b(1:3);
    illum_2 = b(4:6);
    illum_3 = b(7:9);
    
    illum_1 = illum_1/norm(illum_1);
    illum_2 = illum_2/norm(illum_2);
    illum_3 = illum_3/norm(illum_3);
    
    im1_new = [];
    for kk=1:3
        im1_new(kk, :) = (alpha'*E(:,:,kk)*illum_1)'.*corr1;
    end
    im1_new = reshape(im1_new', size(im_nf));

    im2_new = [];
    for kk=1:3
        im2_new(kk, :) = (alpha'*E(:,:,kk)*illum_2)'.*corr2;
    end
    im2_new = reshape(im2_new', size(im_nf));

    im3_new = [];
    for kk=1:3
        im3_new(kk, :) = (alpha'*E(:,:,kk)*illum_3)'.*corr3;
    end
    im3_new = reshape(im3_new', size(im_nf));
    
    results.im1 = im1_new;
    results.im2 = im2_new;
    results.im3 = im3_new;
    
    results.illum1 = illum_1;
    results.illum2 = illum_2;
    results.illum3 = illum_3;
    results.coeff = coeff;
end


