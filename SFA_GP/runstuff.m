% SCRIPT FOR RUNNING Sparse factor analysis with Gaussian process prior

%% Load data:
load RandTestSaureus.mat

% y is a d x N matrix specifying the N observations of the p-dimensional data
[d N] = size(y);

inds_y = ones(size(y));
inds_y = inds_y > 0;


%% set GP prior
% Rescale predictor space for forming correlation matrix:
x = [1:N]./N;

c = 10;  % Gaussian process length scale parameter.  Smaller c will lead to slower changes in correlations.
l = 1;
r = 1e-5;
K = zeros(N);
for ii=1:N
    for jj=1:N
        dist_ii_jj = abs(x(ii)-x(jj));
        K(ii,jj) = l*exp(-c*(dist_ii_jj^2));
    end
end
K = K + diag(r*ones(1,N));
invK = inv(K);
logdetK = 2*sum(log(diag(chol(K))));

prior_params.K.c_prior = 1;
prior_params.K.invK = invK;
prior_params.K.K = K;
prior_params.K.logdetK = logdetK;

%% set noise prior
prior_params.sig.a_sig = 1;
prior_params.sig.b_sig = 0.1;

%% set coefficient prior
prior_params.hypers.a_phi = 1.5;
prior_params.hypers.b_phi = 1.5;
prior_params.hypers.a1 = 2; 
prior_params.hypers.a2 = 2; 

%% settings
settings.t = 10;    % truncation level for basis of latent GPs
settings.k = 10;    % latent factor dimension
settings.Niter = 100;  % number of Gibbs iterations to run
settings.saveEvery = 100;  % how often to write stats to disk
settings.storeEvery = 10;  % how often to store samples in struct
settings.saveMin = 1;
settings.saveDir = ‘../SFA_GP/test';  % name of directory that will be created to save stats in
settings.trial = 1;  % trial number 
settings.sample_K_flag = 3; % choose to fix bandwidth parameter instead of sampling it
settings.latent_mean = 1; % indicate whether one wants to model a non-zero latent mean \mu(x)
settings.inds_y = inds_y;

% Indicate whether we are restarting from a previous chain:
restart = 0; 
if restart
    settings.lastIter = 10000; % last saved iteration we want to load
end

%% call the model
SFA_GP(y,prior_params,settings,restart,true_params);


