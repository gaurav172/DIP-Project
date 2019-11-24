function [qq_density, density]= calDensityMap(qq, numBins)



% qx = (qq(1, :) +1)./2;
qx = qq(1, :);
qx = qx+1;
qx = qx ./ 2;

% qy = (qq(2, :) +1)./2;
qy = qq(2, :);
qy = qy+1;
qy = qy ./ 2;

% qx = (qq(1, :) - min(qq(1, :)))./(max(qq(1, :)) - min(qq(1, :)));
% qy = (qq(2, :) - min(qq(2, :)))./(max(qq(2, :)) - min(qq(2, :)));

xgrid = qx*numBins;
xgrid = xgrid - qx;
xgrid = round(1+qx);
% xgrid = round(1 + qx * (numBins - 1));
% ygrid = round(1 + qy * (numBins - 1));
ygrid = qy*numBins;
ygrid = ygrid - qy;
ygrid = round(1+qy);

idx = sub2ind([numBins numBins], ygrid, xgrid);

density = zeros(numBins, numBins);
idU = unique(idx);
qq_density = zeros(1, size(qq, 2));

for kk = 1:length(idU)
    tempIDD = find(idx == idU(kk));
    density(idU(kk)) = length(tempIDD);
    qq_density(:, tempIDD) = length(tempIDD);
end
