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

%% This function is to implement the RANSAC function


function [n0, nx] = ransac_2d_subspace(pts, inlier_angle, iter, inliers, M)
    
    t1 = size(pts,2)
    
    if ~exist('inlier_angle', 'var')
        inlier_angle = (1*pi/180);
    end
    
    if ~exist('iter', 'var')
        iter = 100;
    end
    
    if ~exist('inliers', 'var')
        inliers = floor(t1/2);
    end


    if ~exist('M', 'var')
        M = 2;
    end
    
    max_inlier = 0;
    g1 = pts(:,x0(1));
    g2 = pts(:,x0(2));



    for kk=1:iter
        x0 = randperm(t1, 2);
        n0 = cross(g1,g2)/norm(cross(g1,g2));

        a1 = n0' * pts;
        a2 = asin(a1);
        ang = abs(a2);

        in_count = length(find(ang <= inlier_angle));
        if (max_inlier <= in_count)
            max_inlier = in_count;
            n0_max = n0;

        end

    end

    a1 = n0' * pts;
    a2 = asin(a1);
    ang = abs(a2);
    % ang = abs(asin(n0_max'*pts));
    % ang = abs(acos(n0_max'*pts));
    in_idx = find(ang <= inlier_angle);

    fprintf('Found %1.3f % \n', max_inlier/t1);

    [U, S, V] = svds(pts(:, in_idx), 3);
    n0 = U(:, 3);
    nx = U(:, 1:2);
