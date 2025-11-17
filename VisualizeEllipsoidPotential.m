function VisualizeEllipsoidPotential
clc, clear, close all

m = 1481.90e20;
a = 2634.40e3;
b = 2634.10e3;
c = 2633.80e3;
G = 6.673e-11;

theta_ = linspace(-pi/2, pi/2, 150);
phi_   = linspace(0, 2*pi, 150);
[theta, phi] = meshgrid(theta_, phi_);
[x, y, z, r] = get_xyzr(theta, phi, a, b, c);
Phi = G*m./r .* (1 + ((2*a^2-b^2-c^2)*x.^2 + (2*b^2-a^2-c^2)*y.^2 + (2*c^2-a^2-b^2)*z.^2) / 10 ./r.^4);
C = contourc(theta_, phi_, Phi, 30);

pos = 1;
N = size(C,2);

figure
hold on
axis equal
view(3)

surf(x, y, z, Phi, 'EdgeColor', 'none')
colormap('jet')
colorbar

while true
  n = C(2, pos);
  theta = C(1, pos+1 : pos+n);
  phi   = C(2, pos+1 : pos+n);
  [x, y, z, r] = get_xyzr(theta, phi, a, b, c);
  plot3(x, y, z, 'k')
  pos = pos + n + 1;
  if pos > N
    break
  end
end

function [x, y, z, r] = get_xyzr(theta, phi, a, b, c)

r = sqrt(1./(cos(theta).^2 .* cos(phi).^2/a^2+cos(theta).^2 .* sin(phi).^2/b^2+sin(theta).^2/c^2));
x = r.*cos(theta).*cos(phi); 
y = r.*cos(theta).*sin(phi); 
z = r.*sin(theta);
