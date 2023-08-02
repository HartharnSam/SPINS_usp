function [error_map, u, coeff] = eof_error_map(data,num_modes)
% EOF_ERROR_MAP computes Empirical Orthogonal Function reconstructions up to and
% including the number of modes specified in num_modes using MATLAB's built in svds.
% The infinity norm error of the reconstruction at each time is stored, for each time,
% and for each reconstruction.
%
% Data needs to be snapshots in columns, so a matrix with M rows giving spatial
% locations and N columns giving timesteps, with time increasing with the column index.
%
%  [error_map] = EOF_ERROR_MAP(data,num_modes)
%
% returns:
% -error_map the time slice error of the eof reconstruction with k modes, where k ranges from
% 1 to num_modes.  Size [num_modes, N] with the k mode reconstruction corresponding to the
% kth row of error_map, and the column index corresponding to time.


%% Setup
[M,N] = size(data); % find dimensions of data
data = bsxfun(@minus, data, mean(data,2)); % remove mean

%% Compute EOFs by svds
[u,~,~]=svds(data/sqrt(N-1),num_modes); % perform the SVD, v is not used.
% the columns of u are the eigenvectors.
% NOTE: data is not transposed as in the Kutz 2013 code on page 394 as this is incorrect, even in his notation

%% Find timeseries coefficients
coeff=u'*data; % project the mean centered data onto the basis, size M x N

%% Calculate Error Map
error_map = zeros(num_modes,N); % error map pre-allocation
recon_temp = zeros(M,N); % loop variable for reconstruction

%% infinity norm Error Map

for ii=1:num_modes % reconstruction to maximum of num_modes
    recon_temp=recon_temp+bsxfun(@times,u(:,ii),coeff(ii,:)); % update reconstruction
    for jj = 1:N % loop through timesteps for current reconstruction
        error_map(ii,jj) = norm(recon_temp(:,jj)-data(:,jj),2); % infinity norm
    end
end
end



