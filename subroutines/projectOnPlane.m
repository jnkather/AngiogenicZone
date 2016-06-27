% created by JN Kather and CA Weis 2015-2016
% jakob.kather@medma.uni-heidelberg.de
% license: MIT license, see separate file for license and disclaimer
%
% parts of the code are based on the following publication
% Kather, JN et al. Continuous representation of tumor microvessel 
% density and detection of angiogenic hotspots in histological 
% whole-slide images. Oncotarget 5, (2015). DOI: 10.18632/oncotarget.4383
%
% others are based on this publication:
% Kather, JN et al. New Colors for Histology: Optimized Bivariate 
% Color Maps Increase Perceptual Contrast in Histological Images. 
% PLoS One 10, e0145572 (2015). DOI: 10.1371/journal.pone.0145572
%
% project points on a plane in 3D
% N=normal, point=Point on Plane
function pointsOut = projectOnPlane(pointsIn, normal, point) 
    pointsIn=double(pointsIn);
    [m,n]=size(pointsIn);
    normal = normal/norm(normal);  % normalize the normal vector
    N2 = normal.'*normal;
    pointsOut = pointsIn*(eye(3)-N2)+repmat(point*N2,m,1);
end