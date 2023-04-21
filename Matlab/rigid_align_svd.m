% Reference rigid registration via the SVD: https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
% p = 3xN real matrix of fiducial positions on MOVING model (ordered/paired
% with FIXED model fiducials).
% q = 3xN real matrix of fiducial positions on FIXED model

function [p_new,TF,RMSE] = rigid_align_svd(p,q)

% center moving points
p_bar = mean(p,2);
X = p - p_bar;

% center fixed points
q_bar = mean(q,2);
Y = q - q_bar;

% compute covariance matrix
W = eye(size(Y,2)); % identity weighting matrix
S = X*W*Y';

% compute SVD
[U,~,V] = svd(S);

% compute homogeneous transformation matrix
W2 = eye(size(U));
W2(end,end) = det(V*U');
R = V*W2*U';
t = q_bar - R*p_bar;
TF = [R, t; zeros(1,3), 1];

% transform and compute residual RMSE (ie. fiducial registration error)
p_new = hTF(p,TF,0);
mse = mean(vecnorm(q-p_new).^2);
RMSE = sqrt(mse); % fiducial registration error

end