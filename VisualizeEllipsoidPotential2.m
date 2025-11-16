function VisualizeEllipsoidPotential2
clc; clear; close all;

%% --- Physical parameters ---
m = 1481.90e20;
a = 2634.40e3;
b = 2634.10e3;
c = 2633.80e3;
G = 6.673e-11;

%% --- Angular grid on ellipsoid ---
theta_ = linspace(-pi/2, pi/2, 150);
phi_   = linspace(0, 2*pi, 150);
[theta, phi] = meshgrid(theta_, phi_);

%% --- Convert (theta,phi) → (x,y,z) on ellipsoid ---
[x, y, z, r] = get_xyzr(theta, phi, a, b, c);

%% --- Gravitational potential (correct formula) ---
Phi = G*m ./ r .* (1 + ...
    ((2*a^2 - b^2 - c^2).*x.^2 + ...
     (2*b^2 - a^2 - c^2).*y.^2 + ...
     (2*c^2 - a^2 - b^2).*z.^2) ./ (10 .* r.^4));

%% --- Contours in parameter space ---
C = contourc(theta_, phi_, Phi, 30);

%% --- PLOT ---
figure('Color','w'); hold on;
axis equal; view(3);

% Surface
surf(x, y, z, Phi, 'EdgeColor','none', 'FaceAlpha', 0.95);
colormap(jet);
shading interp;
colorbar;
title('Ellipsoidal Gravitational Potential','FontSize',14);

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Lighting for better 3D visibility
camlight('headlight');
lighting gouraud;

%% --- Draw contour lines on ellipsoid ---
pos = 1;
total = size(C,2);

while pos < total
    n = C(2,pos);
    if n < 1
        break
    end
    
    theta_c = C(1, pos+1 : pos+n);
    phi_c   = C(2, pos+1 : pos+n);

    [xc,yc,zc,~] = get_xyzr(theta_c, phi_c, a, b, c);
    plot3(xc, yc, zc, 'k', 'LineWidth', 1.0);

    pos = pos + n + 1;
end

end

%% ========================================================================
%                           Coordinate transform
% ========================================================================
function [x, y, z, r] = get_xyzr(theta, phi, a, b, c)
% Planetocentric mapping onto a triaxial ellipsoid

% radius of ellipsoid in given direction
r = sqrt(1 ./ ( ...
     (cos(theta).^2 .* cos(phi).^2)/a^2 + ...
     (cos(theta).^2 .* sin(phi).^2)/b^2 + ...
     (sin(theta).^2)               /c^2 ));

% Cartesian coordinates
x = r .* cos(theta).*cos(phi);
y = r .* cos(theta).*sin(phi);
z = r .* sin(theta);
end
