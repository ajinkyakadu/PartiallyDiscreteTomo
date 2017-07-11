function [m,bgm,mS,bgmS] = createPhantom(x,z,options)
%createPhantom Generates a random phantom model of discrete anomaly in
%       variational background
% Summary:
%
% Input:
%   x is a range vector (horizontal direction)
%   z is a depth vector (vertical direction)
%
% Output:
%   m is a final model containing anomaly and background (a vector 
%           of size n = prod(nx,nz)
%   bgm is the background model of same size as m
%
% Supported Input options:
%   m1 - anomaly parameter value in (default : 1)
%   xwidth - width percentage of the salt in range/x-direction, allowed values 0 to 1, 1 means complete range (default : 0.5)
%   zwidth - width percentage of the salt in depth/z-direction, allowed values 0 to 1, 1 means complete depth (default : 0.5)
%   xoffset - offset of salt in range/x-direction in m, can take +ve and -ve values, maximum is half of range (default : 0)
%   zoffset - offset of salt in depth/z-direction in m, can take +ve and -ve values, maximum is half of depth (default : 0)
%   nrand - number of points in the valid domain to generate the boundary of salt (default : 20)
%   randseed - random seed for generating radom points to form boundary (default : 0)
%
%
% Author: Ajinkya Kadu
%         Mathematical Institute,
%         Utrecht University, The Netherlands
%
% Date : January 2017

if nargin < 3
    options = [];
end

m1      = getoptions(options,'m1',1);
xwidth  = getoptions(options,'xwidth',0.5);
zwidth  = getoptions(options,'zwidth',0.5);
xoffset = getoptions(options,'xoffset',0);
zoffset = getoptions(options,'zoffset',0);
nrand   = getoptions(options,'nrand',10);
ri      = getoptions(options,'randseed',0);
smooth  = getoptions(options.bg,'smooth',10);
bmax    = getoptions(options.bg,'bmax',0.5);
type    = getoptions(options,'type',1);
rseed   = getoptions(options,'rseed',1);
gV      = getoptions(options,'gV',[]);

nz = length(z);
nx = length(x);

[zz,xx] = ndgrid(z,x);

s = RandStream('mt19937ar','Seed',rseed);
RandStream.setGlobalStream(s);


m0 = 1.5*rand(nz,nx);
% m0 = vsmooth(m0,bgsmooth);
H = fspecial('disk',smooth);
m0 = imfilter(m0,H,'replicate');
m0 = m0 - min(m0(:));
m0 = (bmax/max(m0(:)))*m0;

bgm = m0(:);

% segmented background for DART
bgmS = bgm(:);
gvl = gV(1);
gvr = (gV(1)+gV(2))/2;
bgmS((bgm >= gvl) & (bgm <gvr)) = gV(1);
for i=2:(length(gV)-1)
    gvl = (gV(i-1)+gV(i))/2;
    gvr = (gV(i)+gV(i+1))/2;
    bgmS((bgm >= gvl) & (bgm <gvr) ) = gV(i);
end
gvl = (gV(end-1)+gV(end))/2;
gvr = gV(end);
bgmS((bgm >= gvl) & (bgm <gvr)) = gV(end);

bgm = reshape(bgm,nz,nx);
bgmS= reshape(bgmS,nz,nx);

m = m0;
mS = bgmS;

switch type
    case 1
        rng(ri,'twister');
        
        % random points spread over the computational domain
        xt = (min(x)+max(x))/2 + xwidth*(max(x)-min(x))*(rand(nrand,1)-0.5) + xoffset;
        zt = (min(z)+max(z))/2 + zwidth*(max(z)-min(z))*(rand(nrand,1)-0.5) + zoffset;
        
        % get points on the boundary
        k = boundary(xt,zt);
        xv = xt(k);
        zv = zt(k);
        
        % locate points inside the boundary and set them to salt velocity
        [in,~] = inpolygon(xx,zz,xv,zv);
        m(in) = m1;
        mS(in)= m1;
    
    case 2
        I = imread([pwd '/phantoms/foam_128.png']);
        n = size(I,1);
        I = double(I(:)/max(I(:)));
        m(I == 1) = m1;
        mS(I==1)  = m1;
        
    case 3
        I = imread([pwd '/phantoms/molecule_128.png']);
        n = size(I,1);
        I = double(I(:)/max(I(:)));
        m(I==1) = m1;
        mS(I==1)= m1;
end

end