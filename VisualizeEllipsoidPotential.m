function VisualizeEllipsoidPotential
clc, clear, close all

%% Physical parameters
m = 1481.90e20;         % Mass of Ganymede (kg)
a = 2634.40e3;          % Semi-axis along x (m)
b = 2634.10e3;          % Semi-axis along y (m)
c = 2633.80e3;          % Semi-axis along z (m)
G = 6.673e-11;          % Gravitational constant (m³/kg/s²)
omega = 2*pi/(7.15*24*3600);  % Angular velocity of Ganymede (rad/s), period ~7.15 days

%% Angular grid
theta_ = linspace(-pi/2, pi/2, 150);
phi_   = linspace(0, 2*pi, 150);
[theta, phi] = meshgrid(theta_, phi_);

%% Convert to Cartesian coordinates on ellipsoid
[x, y, z, r] = get_xyzr(theta, phi, a, b, c);

%% Calculate potentials
% Gravitational potential (quadrupole approximation)
Phi_grav = -G*m./r .* (1 + ((2*a^2-b^2-c^2)*x.^2 + (2*b^2-a^2-c^2)*y.^2 + (2*c^2-a^2-b^2)*z.^2) / 10 ./r.^4);

% Centrifugal potential
Phi_centrifugal = -0.5 * omega^2 * (x.^2 + y.^2);

% Effective potential
Phi_eff = Phi_grav + Phi_centrifugal;

%% Contour lines for gravitational potential
C = contourc(theta_, phi_, Phi_grav, 30);

%% Plotting
pos = 1;
N = size(C,2);

figure('Color','w')
hold on
axis equal
view(3)

% Gray ellipsoid base
surf(x, y, z, ones(size(z)), 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3)

% Effective potential colored surface
scale_factor = 1.002;
surf(scale_factor*x, scale_factor*y, scale_factor*z, Phi_eff, 'EdgeColor', 'none', 'FaceAlpha', 0.7)

% Blue to yellow colormap
colormap([linspace(0,1,64)', linspace(0,1,64)', linspace(1,0,64)'])
colorbar
title('Ganymede: Gravitational Field Model with Rotation')
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

%% Draw dashed contour lines
while true
  n = C(2, pos);
  theta = C(1, pos+1 : pos+n);
  phi   = C(2, pos+1 : pos+n);
  [x, y, z, r] = get_xyzr(theta, phi, a, b, c);
  plot3(x, y, z, 'k--', 'LineWidth', 1)
  pos = pos + n + 1;
  if pos > N
    break
  end
end

%% Add gravitational force arrows
theta_arrow = linspace(-pi/2, pi/2, 15);
phi_arrow = linspace(0, 2*pi, 20);
[theta_a, phi_a] = meshgrid(theta_arrow, phi_arrow);
[x_a, y_a, z_a, r_a] = get_xyzr(theta_a, phi_a, a, b, c);

% Gravitational force direction (centripetal, pointing toward center)
arrow_scale = max(a,max(b,c)) * 0.05;
base_offset = 1.02;

for i = 1:numel(x_a)
    x_base = x_a(i) * base_offset;
    y_base = y_a(i) * base_offset;
    z_base = z_a(i) * base_offset;
    
    quiver3(x_base, y_base, z_base, ...
            x_a(i) - x_base, y_a(i) - y_base, z_a(i) - z_base, ...
            0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end

%% Helper function
function [x, y, z, r] = get_xyzr(theta, phi, a, b, c)

r = sqrt(1./(cos(theta).^2 .* cos(phi).^2/a^2+cos(theta).^2 .* sin(phi).^2/b^2+sin(theta).^2/c^2));
x = r.*cos(theta).*cos(phi); 
y = r.*cos(theta).*sin(phi); 
z = r.*sin(theta);