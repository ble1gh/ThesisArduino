FV.faces = [5 3 1; 3 2 1; 3 4 2; 4 6 2];
FV.vertices = [2.5 8.0 1; 6.5 8.0 2; 2.5 5.0 1; 6.5 5.0 0; 1.0 6.5 1; 8.0 6.5 1.5];
points = [2 4 2; 4 6 2; 5 6 2];
[distances,surface_points] = point2trimesh(FV, 'QueryPoints', points);
patch(FV,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
%%
plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
plot3M(points,'*r')
%plot3M(surface_points,'*k')
%plot3M(reshape([shiftdim(points,-1);shiftdim(surface_points,-1);shiftdim(points,-1)*NaN],[],3),'k')
