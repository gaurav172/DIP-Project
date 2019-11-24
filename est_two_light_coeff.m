%% This function is to implement the esimtation of lighting coefficients

% Input
% 1. gamma: 3 by  M*N
% 2. mask_t: mask for the image M by N
% 3. cfactor: the threshold to find the corner points

%   

% Output
% 1. illum1, illum2: estimated lighting coefficents

function [illum1, illum2] = est_two_light_coeff(gamma, mask_t, cfactor)
    [~, nx] = ransac_2d_subspace(gamma(:, mask_t>0), (pi*0.3)/180, 1000);
    theta = ((0:.1:360)*pi)/180;
    len_t = length(theta)
    Circle = nx*[cos(theta); sin(theta)];

    iit = find(mask_t > 0);
    for ii = 1:length(iit)
        Circle2 = Circle'*gamma(:, iit(ii));
        [~, idx(ii)] = max(Circle2);
    end
    [h_idx] = hist(idx, 1:size(Circle, 2));

    cutoff = cfactor*sum(h_idx)/len_t;

    tmp = find(h_idx > cutoff);
    theta1 = theta(min(tmp));
    theta2 = theta(max(tmp));
    illum1 = nx*[cos(theta1); sin(theta1)];
    illum2 = nx*[cos(theta2); sin(theta2)];

