clf; hold on;
N=100;
C=zeros(N,N,3); C(:,:,1) = 0; C(:,:,2) = 1; C(:,:,3) = 0;
[X,Y,Z] = sphere(N);
surf(X,Y,Z,C,'EdgeAlpha',0,'FaceLighting','phong');
X = 0+0.7*(X-0); Y = 2+0.7*(Y-2); Z = 2+0.7*(Z-2);
surf(X,Y,Z,C,'EdgeAlpha',0,'FaceLighting','phong');
axis equal; axis off; camlight right;
hold off;
