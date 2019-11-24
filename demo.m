% set up the route
addpath(genpath('extern'))
addpath('utilits/')

% load the data to process
disp('Loading data...')
load('data/reflectance_illum_camera.mat')
load('data/images/two_lights/book.mat')
disp('Done.')

% generate basis for both reflectance/illumination
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Generate basis....')
[UR, UL] = genBase(L, C, R, 'wpca1');
for ind = 1:3
    E(:,:,ind) = UR'*diag(C(:, ind))*UL;
end
opt.E = E;

flash_light = .025*ones(size(L(:, 1)));
f = UL'*flash_light;
f = f/norm(f);
opt.f = f;

%% Demo for TWO lights

%% These parameter can be adjusted for better performance
opt.lambda = 1e-5;
opt.cutoff = .1;
opt.color_correct = 1;
opt.light_number = 2;
opt.shadow_mask = 1 - mask;
disp('Done.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
disp('Processing the data....')
results = solve_light_sep(im_nf, im_f, mask, opt);
disp('Done');

%% Visulize the results
subplot(2,1,1);
imshow([im_nf im_f].^(1/2.2))
title('Image - Flash and No-flash');
subplot(2,1,2);
imshow([results.im1 results.im2].^(1/2.2));
title('Seperated images');
fprintf('Program paused. Press enter to continue.\n');
pause;

%% Demo for Three lights

%% These parameter can be adjusted for better performance
disp('Loading data...')
load('data/images/three_lights/two_paintings.mat')
disp('Done.')

opt.lambda = 1e-5;
opt.cutoff = .1;
opt.color_correct = 1;
opt.light_number = 3;
opt.shadow_mask = shadow_mask;
disp('Done.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
disp('Processing the data....')
results = solve_light_sep(im_nf, im_f, mask, opt);
disp('Done');

%% Visulize the results
subplot(2,1,1);
imshow([im_nf im_f].^(1/2.2))
title('Image - Flash and No-flash');
subplot(2,1,2);
imshow([results.im1 results.im2 results.im3].^(1/2.2));
title('Seperated images');