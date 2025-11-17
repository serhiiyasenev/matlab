function VisualizeEllipsoidPotential2
clc; clear; close all;

%% --- Physical parameters ---
m = 1481.90e20;         % Mass of Ganymede (kg)
a = 2634.40e3;          % Semi-axis along x (m)
b = 2634.10e3;          % Semi-axis along y (m)
c = 2633.80e3;          % Semi-axis along z (m)
G = 6.673e-11;          % Gravitational constant (m³/kg/s²)
omega = 2*pi/(7.15*24*3600);  % Angular velocity of Ganymede (rad/s), period ~7.15 days

%% --- Angular grid on ellipsoid ---
theta_ = linspace(-pi/2, pi/2, 150);
phi_   = linspace(0, 2*pi, 150);
[theta, phi] = meshgrid(theta_, phi_);

%% --- Convert (theta,phi) → (x,y,z) on ellipsoid ---
[x, y, z, r] = get_xyzr(theta, phi, a, b, c);

%% --- Gravitational potential (quadrupole approximation) ---
Phi_grav = -G*m ./ r .* (1 + ...
    ((2*a^2 - b^2 - c^2).*x.^2 + ...
     (2*b^2 - a^2 - c^2).*y.^2 + ...
     (2*c^2 - a^2 - b^2).*z.^2) ./ (10 .* r.^4));

%% --- Centrifugal potential (due to rotation) ---
% U_centrifugal = -0.5 * omega^2 * rho^2, where rho^2 = x^2 + y^2
Phi_centrifugal = -0.5 * omega^2 * (x.^2 + y.^2);

%% --- Effective potential (gravitational + centrifugal) ---
Phi_eff = Phi_grav + Phi_centrifugal;

%% --- Contours in parameter space for gravitational potential ---
C = contourc(theta_, phi_, Phi_grav, 30);

%% --- PLOT ---
figure('Color','w'); hold on;
axis equal; view(3);

% Gray ellipsoid surface
surf(x, y, z, ones(size(z)), 'EdgeColor','none', ...
     'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3);

% Effective potential overlay (colored surface slightly above gray surface)
scale_factor = 1.002;  % Slightly offset to avoid z-fighting
surf(scale_factor*x, scale_factor*y, scale_factor*z, Phi_eff, ...
     'EdgeColor','none', 'FaceAlpha', 0.7);

% Blue to yellow colormap for effective potential
colormap([linspace(0,1,64)', linspace(0,1,64)', linspace(1,0,64)']);
shading interp;
cb = colorbar;
ylabel(cb, 'Effective Potential [J/kg]', 'FontSize', 10);
title('Ganymede: Gravitational Field Model with Rotation', 'FontSize', 14);

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Lighting for better 3D visibility
camlight('headlight');
lighting gouraud;

%% --- Draw contour lines on ellipsoid (dashed black lines) ---
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
    plot3(xc, yc, zc, 'k--', 'LineWidth', 1.0);  % Dashed black lines

    pos = pos + n + 1;
end

%% --- Add gravitational force vectors (red arrows) ---
% Sample points for arrows (reduced grid for visibility)
theta_arrow = linspace(-pi/2, pi/2, 15);
phi_arrow   = linspace(0, 2*pi, 20);
[theta_a, phi_a] = meshgrid(theta_arrow, phi_arrow);

[x_a, y_a, z_a, r_a] = get_xyzr(theta_a, phi_a, a, b, c);

% Calculate gravitational force direction (pointing toward center)
% This is the unit vector from surface point to center
Fx = -x_a ./ r_a;
Fy = -y_a ./ r_a;
Fz = -z_a ./ r_a;

% Normalize and scale for visualization
arrow_scale = max(a,max(b,c)) * 0.05;  % Arrow length
Fx_norm = Fx * arrow_scale;
Fy_norm = Fy * arrow_scale;
Fz_norm = Fz * arrow_scale;

% Draw arrows: base slightly above surface, tip on surface
base_offset = 1.02;
for i = 1:numel(x_a)
    x_base = x_a(i) * base_offset;
    y_base = y_a(i) * base_offset;
    z_base = z_a(i) * base_offset;
    
    x_tip = x_a(i);
    y_tip = y_a(i);
    z_tip = z_a(i);
    
    % Draw arrow from base to tip
    quiver3(x_base, y_base, z_base, ...
            x_tip - x_base, y_tip - y_base, z_tip - z_base, ...
            0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
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
