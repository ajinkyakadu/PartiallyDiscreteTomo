% addpath(genpath('/path/to/astra-toolbox/astra-1.8/'));
% addpath(genpath('/path/to/spot/operator/spot/'));

fprintf('Please add the paths to ASTRA toolbox and SPOT operator in startup.m \n');
fprintf('ASTRA toolbox available here: https://github.com/astra-toolbox/astra-toolbox \n');
fprintf('SPOT is available here: https://github.com/mpf/spot \n');

addpath(genpath([pwd '/mbin/']));
addpath(genpath([pwd '/phantoms/']));

fprintf('If ASTRA and SPOT added, then success! \n');