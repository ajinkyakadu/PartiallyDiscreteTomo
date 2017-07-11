function [h,d] = heavi(x,options)
%
% heavi Heaviside function (approximated) operating on vector x
%
% Summary
%
% Input:
%   x   - Real-valued vector of any size
%   options - is a struct containing parameters (defaults are used for non-existent or blank fields)
%
% Output:
%   h   = vector of heaviside values of size of x
%   d   = vector of delta (derivative of heaviside) values of size of x
% 
% Supported Input options:
%   type - Type of Heaviside [ global | (compact) ]
%   epsi - Heaviside epsilon (default : 0.01)
%   thr  - Threshold for level-set, i.e x = 0 (default : 0)


if nargin<2
    options = [];
end

type    = getoptions(options,'type','compact_new');
epsi    = getoptions(options,'epsi',0.1);
thr     = getoptions(options,'thr',0);
k       = getoptions(options,'k',0);

% x = x - (min(x)+max(x))/2;

x = x - thr;

switch type
    case 'global'
        h = 0.5*(1 + 2/pi*atan(pi*x/epsi));
        d = 1./(epsi*((x.^2*pi^2)/epsi^2 + 1));
        
    case 'compact'
        h = 0*x;
        d = 0*x;
        id = find((x < epsi) & (x > -epsi));
        h(id) = 0.5*(1 + x(id)/epsi + 1/pi*sin(pi*x(id)/epsi));
        h(x >= epsi) = 1;
        h(x <=-epsi) = 0;
        d(id) = 0.5*(1/epsi)*(1 + cos(pi*x(id)/epsi));
        
    case 'compact_new'
        h = 0*x;
        d = 0*x;
        L = epsi;
        
        id1 = find((x>-L) & (x <= (2*k-1)*L));
        d(id1) = 0.25/((1-k)*L) * (1 + 1/(k*L)*(x(id1)+(1-k)*L) + 1/pi*sin(pi*(x(id1)+(1-k)*L)/(k*L)));
        h(id1) = 0.25/((1-k)*L) * ( x(id1) + 0.5/(k*L)*(x(id1) + (1-k)*L).^2 -(k*L)/pi^2*(cos(pi*(x(id1)+(1-k)*L)/(k*L)))...
            + L - 0.5*(k*L) - (k*L)/pi^2);
        
        id2 = find((x<(1-2*k)*L) & (x>(2*k-1)*L));
        d(id2) = 0.5/((1-k)*L);
        h(id2) = 0.5*(x(id2)+(1-2*k)*L)/((1-k)*L) + 0.5*k/(1-k);
        
        id3 = find((x>=(1-2*k)*L) & (x<L));
        d(id3) = 0.25*(1 - (x(id3)-(1-k)*L)/(k*L) - 1/pi*sin(pi*(x(id3)-(1-k)*L)/(k*L)))/((1-k)*L);
        h(id3) = 0.25/((1-k)*L)*(x(id3) - 0.5/(k*L)*(x(id3)-(1-k)*L).^2 + (k*L)/pi^2*cos(pi*(x(id3)-(1-k)*L)/...
            (k*L)) - (1-2*k)*L + 0.5*(k*L) + (k*L)/pi^2) + 0.5*(2-3*k)/(1-k);
        
        h(x >= L) = 1;
        
end

end