% apply a homogeneous transformation matrix to a set of points
% this function simply reduces the need to augment the point matrix
% with a column of ones and strip it away after applying the transformation
%
% inPoints:     m x n matrix of points to transform
% T:            homogeneous transformation matrix (typically 3x3 for 2D or 4x4 for 3D)
% premult_flag: (optional) 0: postmultiply T by points; 1: premultiply T by points
function outPoints = hTF( inPoints, T, premult_flag )

% determine whether we are set up to pre- or post-multiply the point
% matrix by the transformation matrix T
lastColZeros = prod(T(1:end-1,end) == zeros(size(T,1)-1,1));
lastRowZeros = prod(T(end,1:end-1) == zeros(1,size(T,2)-1));
if( nargin < 2 )
    error('Too few input arguments!');
elseif( nargin == 2)
    if( lastColZeros && ~lastRowZeros )
        premult_flag = 1;
    elseif( lastRowZeros && ~lastColZeros )
        premult_flag = 0;
    else
        error('Cannot determine whether to pre- or post- multiply, please specify explicitly with premult_flag argument!');
    end
end

% perform transformation
% TODO: this assumes that inPoints are given correctly as row or column vectors
switch(premult_flag)
    case 0
        % post-multiply matrix by points to transform
        outPoints = T*[inPoints; ones(1,size(inPoints,2))];
        outPoints = outPoints(1:end-1,:);
    case 1
        % pre-multiply matrix by points to transform
        outPoints = [inPoints ones(size(inPoints,1),1)]*T;
        outPoints = outPoints(:,1:end-1);
    otherwise
        error('Invalid value for premult_flag');
end