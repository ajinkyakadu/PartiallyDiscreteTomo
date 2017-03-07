function show(img, varargin)

warning('off', 'Images:initSize:adjustingMag');

% figure handle
h = gcf;
if h < 1
    error('Could not obtain figure handle!');
end

figure(h);

isDocked = strcmp('docked', get(h, 'WindowStyle'));

% if ~hasDisplay
%     return
% end

if isvector(img)
    N = numel(img);
    n = sqrt(N);
    if round(n)^2 == N
        img = reshape(img,n,n);
    else
        plot(img, varargin{:});
        return
    end
end

% prevent warning
if isDocked
    varargin{1,end+1} = 'in';
    varargin{1,end+1} = 'fit';
end

imshow(squeeze(img), [], varargin{:});
drawnow;
