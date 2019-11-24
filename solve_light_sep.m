


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
if isfield(opt,'pMask')
    pMask = opt.pMask;
end

lambda = 1e-3;
if isfield(opt,'lambda')
    lambda = opt.lambda;
end

cfactor = 1;
if isfield(opt,'cutoff')
    cfactor = opt.cutoff;
end

E = opt.E;

f = opt.f;

Q = opt.light_number;

siz = size(im_nf);

%% generate difference image
D_img = im_f - im_nf;
D_img(D_img<0) = 0;
Temp_Mat = repmat(mask,[1 1 3]);
diff_img = D_img.*Temp_Mat;

for i=1:3
    Amat(:,i) = E(:,:,i)*f;
end


diff_img_vec = reshape(diff_img, [], 3);
diff_img_vec = diff_img_vec';
img_nf_vec = reshape(im_nf, [], 3)';
img_nf_vec = img_nf_vec;

%% get alpha image
alpha = (Amat');
alpha = alpha\diff_img_vec;

if ~isfield(opt, 'shadow_mask')
    alpha = shadowRM(im_nf, im_f, mask, reshape(alpha', siz), opt.shadow_mask)';
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
        Bmat = Bmat/(1e-10+norm_alpha(kk));
        
        
        nf_intensity = img_nf_vec(:, kk);
        num_1 = Bmat'*Bmat;
        num_2 = lambda*eye(3);
        deno = Bmat'*nf_intensity;
        z = 0;
        if ~pMask(kk)
            z = 5;
        else
            z = 1;
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
s_img(:,:,1) = imguidedfilter(s_img(:,:,1), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img(:,:,2) = imguidedfilter(s_img(:,:,2), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img(:,:,3) = imguidedfilter(s_img(:,:,3), im_nf_t, 'DegreeOfSmoothing', smoothness*diff(getrangefromclass(im_nf)).^2);
s_img = reshape(s_img, [], 3)';
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
    [illum1,illum2] = est_two_light_coeff(gamma, mask_t, cfactor);
    illum = [illum1,illum2]
    L_mat = illum;
    %% reproject the shadings

    for i=1:2
        mR(i,:) = (alpha'*E(:,:,1)*illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end

    for i=1:2
        mG(i,:) = (alpha'*E(:,:,2)*illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end

    for i=1:2
        mB(i,:) = (alpha'*E(:,:,3)*illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end
    for kk = 1:size(img_nf_vec,2)
        Amat_ = [mR(1,kk) mR(2,kk); mG(1,kk) mG(2,kk); mB(1,kk) mB(2,kk)];
        bvec_ = img_nf_vec(:, kk);
        coeff(:, kk) = Amat_\bvec_;
    end

%%
    for i=1:2
        corr(i,:) = coeff(i,:).*beta_norm./(1e-10+norm_alpha);
    end

    dc1 = repmat(corr(1,:), [3,1]);
    dc2 = repmat(corr(2,:), [3,1]);   
    nc1 = (alpha'*E(:,:,1))';
    nc2 = (alpha'*E(:,:,2))';
    nc3 = (alpha'*E(:,:,3))';

    MR1 = nc1.*dc1;
    MR2 = nc1.*dc2;

    MG1 = nc2.*dc1;
    MG2 = nc2.*dc2;


    MB1 = nc3.*dc1;
    MB2 = nc3.*dc2;

    b= [MR1' MR2'; MG1' MG2'; MB1' MB2']\(reshape(img_nf_vec', [], 1));

    Illum_n = [b(1:3),b(4:6)];

    Illum_n(:,1) = Illum_n(:,1)/norm(Illum_n(:,1));
    Illum_n(:,2) = Illum_n(:,2)/norm(Illum_n(:,2));
    im1_new = [];
    im1_new(1, :) = (alpha'*E(:,:,1)*Illum_n(:,1))'.*corr(1,:);
    im1_new(2, :) = (alpha'*E(:,:,2)*Illum_n(:,1))'.*corr(1,:);
    im1_new(3, :) = (alpha'*E(:,:,3)*Illum_n(:,1))'.*corr(1,:);

    im1_new = reshape(im1_new', size(im_nf));

    im2_new = [];
    im2_new(1, :) = (alpha'*E(:,:,1)*Illum_n(:,2))'.*corr(2,:);
    im2_new(2, :) = (alpha'*E(:,:,2)*Illum_n(:,2))'.*corr(2,:);
    im2_new(3, :) = (alpha'*E(:,:,3)*Illum_n(:,2))'.*corr(2,:);

    im2_new = reshape(im2_new', size(im_nf));


    im1_wb = [];
    im1_wb(1, :) = (alpha'*E(:,:,1)*f)'.*corr(1,:);
    im1_wb(2, :) = (alpha'*E(:,:,2)*f)'.*corr(1,:);
    im1_wb(3, :) = (alpha'*E(:,:,3)*f)'.*corr(1,:);

    im1_wb = reshape(im1_wb', size(im_nf));

    im2_wb = [];
    im2_wb(1, :) = (alpha'*E(:,:,1)*f)'.*corr(2,:);
    im2_wb(2, :) = (alpha'*E(:,:,2)*f)'.*corr(2,:);
    im2_wb(3, :) = (alpha'*E(:,:,3)*f)'.*corr(2,:);

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
    Illum = [illum1,illum2,illum3];
    
    for i=1:3
        mR(i,:) = (alpha'*E(:,:,1)*Illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end


    for i=1:3
        mG(i,:) = (alpha'*E(:,:,2)*Illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end

    for i=1:3
        mB(i,:) = (alpha'*E(:,:,3)*Illum(:,i))'.*(beta_norm./(1e-10+norm_alpha));
    end

    
    for kk = 1:size(img_nf_vec,2)
        Amat_ = [mR(1,kk) mR(2,kk) mR(3,kk); mG(1,kk) mG(2,kk) mG(3,kk); mB(1,kk) mB(2,kk) mB(3,kk)];
        bvec_ = img_nf_vec(:, kk);
        coeff(:, kk) = Amat_\bvec_;
    end    
        

    for i=1:3
        corr(i,:) = coeff(i,:).*beta_norm./(1e-10+norm_alpha);
    end
    %%
    nc1 = (alpha'*E(:,:,1))';
    nc2 = (alpha'*E(:,:,2))';
    nc3 = (alpha'*E(:,:,3))';
    
    ncc1 = repmat(corr(1,:), [3 1]);
    ncc2 = repmat(corr(2,:), [3 1]);
    ncc3 = repmat(corr(3,:), [3 1]);
    
    MR1 = nc1.*ncc1;
    MR2 = nc1.*ncc2;
    MR3 = nc1.*ncc3;
    
    MG1 = nc2.*ncc1;
    MG2 = nc2.*ncc2;
    MG3 = nc2.*ncc3;
    
    MB1 = nc3.*ncc1;
    MB2 = nc3.*ncc2;
    MB3 = nc3.*ncc3;
    
    b= [MR1' MR2' MR3'; MG1' MG2' MG3'; MB1' MB2' MB3']\(reshape(img_nf_vec', [], 1));

    Illum = [b(1:3),b(4:6),b(7:9)];
    for i=1:3
        Illum(:,i) = Illum(:,i)/norm(Illum(:,i));
    end
    
    im1_new = [];
    im1_new(1, :) = (alpha'*E(:,:,1)*Illum(:,1))'.*corr(1,:);
    im1_new(2, :) = (alpha'*E(:,:,2)*Illum(:,1))'.*corr(1,:);
    im1_new(3, :) = (alpha'*E(:,:,3)*Illum(:,1))'.*corr(1,:);

    im1_new = reshape(im1_new', size(im_nf));

    im2_new = [];
    im2_new(1, :) = (alpha'*E(:,:,1)*Illum(:,2))'.*corr(2,:);
    im2_new(2, :) = (alpha'*E(:,:,2)*Illum(:,2))'.*corr(2,:);
    im2_new(3, :) = (alpha'*E(:,:,3)*Illum(:,2))'.*corr(2,:);

    im2_new = reshape(im2_new', size(im_nf));

    im3_new = [];
    im3_new(1, :) = (alpha'*E(:,:,1)*Illum(:,3))'.*corr(3,:);
    im3_new(2, :) = (alpha'*E(:,:,2)*Illum(:,3))'.*corr(3,:);
    im3_new(3, :) = (alpha'*E(:,:,3)*Illum(:,3))'.*corr(3,:);

    im3_new = reshape(im3_new', size(im_nf));
    
    results.im1 = im1_new;
    results.im2 = im2_new;
    results.im3 = im3_new;
    
    results.illum1 = Illum(:,1);
    results.illum2 = Illum(:,2);
    results.illum3 = Illum(:,3);
    results.coeff = coeff;
end


