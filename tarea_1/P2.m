%%
% a)
% Patrón en 3D
clc;
clear all;

% Parámetros
theta = linspace(0, pi, 720);  % Variando de 0 a pi (ángulo de elevación)
phi = linspace(0, 2*pi, 1440);  % Variando de 0 a 2*pi (ángulo azimutal)

[Theta, Phi] = meshgrid(theta, phi);  % Crea una malla

% Patrón de radiación
R = sin(Theta);

% Convertir a coordenadas cartesianas
X = R .* cos(Phi) .* sin(Theta);
Y = R .* sin(Phi) .* sin(Theta);
Z = R .* cos(Theta);

% Graficar
figure;
surf(X,Y,Z,R);
title('Patrón de radiación en 3D de sin(\theta)');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% Patrón en 2D

% Graficar
figure;
contour(X,Y,R);
title('Corte 2D del patrón de radiación de sin(\theta)');
xlabel('X');
ylabel('Y');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% HPBW

% Calcular el patrón de radiación
P = sin(theta);

% Encontrar el ángulo en el que la potencia es máxima
[max_intensity, max_index] = max(P);
max_angle = theta(max_index);

% Encontrar los ángulos en los cuales la potencia cae a la mitad del máximo
half_max_intensity = max_intensity / 2;
left_index = find(P(1:max_index) <= half_max_intensity, 1, 'last');
right_index = max_index + find(P(max_index:end) <= half_max_intensity, 1, 'first');

% Calcular el HPBW
hpbw = abs(theta(right_index) - theta(left_index));

disp(['HPBW:', num2str(hpbw)]);

% Graficar el patrón de radiación y las líneas para HPBW
figure;
plot(theta, P, 'b');
hold on;
plot([theta(left_index), theta(left_index)], [0, max_intensity], 'r--');
plot([theta(right_index), theta(right_index)], [0, max_intensity], 'r--');
hold off;
title('Patrón de Radiación sin(\theta)');
xlabel('\theta');
ylabel('Intensidad');
grid on;
legend('Patrón de Radiación', 'Límites HPBW');

% Directividad

% Definir la función a integrar (ángulo sólido)
integrand = @(x, y) sin(x).^3; %sin(x);

% Definir los límites de integración para x (primer integral)
x_lower_limit = 0;
x_upper_limit = pi;

% Definir los límites de integración para y (segunda integral)
y_lower_limit = 0;
y_upper_limit = 2 * pi;

% Calcular la doble integral
result = integral2(integrand, x_lower_limit, x_upper_limit, y_lower_limit, y_upper_limit);

fprintf('Resultado de la doble integral: %.6f\n', result);

fprintf('Máxima Directividad: %.6f\n', 4*pi/result);

%%
% b)
% Patrón en 3D
clc;
clear all;

% Parámetros
theta = linspace(0, pi, 720);  % Variando de 0 a pi (ángulo de elevación)
phi = linspace(0, 2*pi, 1440);  % Variando de 0 a 2*pi (ángulo azimutal)

[Theta, Phi] = meshgrid(theta, phi);  % Crea una malla

% Patrón de radiación
R = sin(Theta).^3;

% Convertir a coordenadas cartesianas
X = R .* cos(Phi) .* sin(Theta);
Y = R .* sin(Phi) .* sin(Theta);
Z = R .* cos(Theta);

% Graficar
figure;
surf(X,Y,Z,R);
title('Patrón de radiación en 3D de sin^3(\theta)');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% Patrón en 2D

% Graficar
figure;
contour(X,Y,R);
title('Corte 2D del patrón de radiación de sin^3(\theta)');
xlabel('X');
ylabel('Y');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% HPBW

% Calcular el patrón de radiación
P = sin(theta).^3;

% Encontrar el ángulo en el que la potencia es máxima
[max_intensity, max_index] = max(P);
max_angle = theta(max_index);

% Encontrar los ángulos en los cuales la potencia cae a la mitad del máximo
half_max_intensity = max_intensity / 2;
left_index = find(P(1:max_index) <= half_max_intensity, 1, 'last');
right_index = max_index + find(P(max_index:end) <= half_max_intensity, 1, 'first');

% Calcular el HPBW
hpbw = abs(theta(right_index) - theta(left_index));

disp(['HPBW:', num2str(hpbw)]);

% Graficar el patrón de radiación y las líneas para HPBW
figure;
plot(theta, P, 'b');
hold on;
plot([theta(left_index), theta(left_index)], [0, max_intensity], 'r--');
plot([theta(right_index), theta(right_index)], [0, max_intensity], 'r--');
hold off;
title('Patrón de Radiación sin^3(\theta)');
xlabel('\theta');
ylabel('Intensidad');
grid on;
legend('Patrón de Radiación', 'Límites HPBW');

% Directividad

% Definir la función a integrar (ángulo sólido)
integrand = @(x, y) sin(x).^7; %sin(x);

% Definir los límites de integración para x (primer integral)
x_lower_limit = 0;
x_upper_limit = pi;

% Definir los límites de integración para y (segunda integral)
y_lower_limit = 0;
y_upper_limit = 2 * pi;

% Calcular la doble integral
result = integral2(integrand, x_lower_limit, x_upper_limit, y_lower_limit, y_upper_limit);

fprintf('Resultado de la doble integral: %.6f\n', result);

fprintf('Máxima Directividad: %.6f\n', 4*pi/result);

%%
% c)
% Patrón en 3D
clc;
clear all;

% Parámetros
alpha = 1.0;
beta = 10.0;
theta = linspace(-pi/2, pi/2, 720);  % Variando de -pi/2 a pi/2 (ángulo de elevación)
phi = linspace(0, 2*pi, 1440);  % Variando de 0 a 2*pi (ángulo azimutal)

[Theta, Phi] = meshgrid(theta, phi);  % Crea una malla

% Patrón de radiación
R = 2.*besselj(1,beta.*alpha.*sin(Theta))./(beta.*alpha.*sin(Theta));

% Convertir a coordenadas cartesianas
X = R .* cos(Phi) .* sin(Theta);
Y = R .* sin(Phi) .* sin(Theta);
Z = R .* cos(Theta);

% Graficar
figure;
surf(X,Y,Z,R);
title('Patrón de radiación en 3D de 2 * J1(βαsin(θ)) / (βαsin(θ))');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% Patrón en 2D

% Graficar

figure;
contour(X,Y,R,50);
title('Corte 2D del patrón de radiación de 2 * J1(βαsin(θ)) / (βαsin(θ))');
xlabel('X');
ylabel('Y');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

figure;
contour(X,Z,R,50);
title('Corte 2D del patrón de radiación de 2 * J1(βαsin(θ)) / (βαsin(θ))');
xlabel('X');
ylabel('Z');
axis equal;  % Para hacer que los ejes tengan la misma escala

shading interp;
colormap('jet');  % Color
colorbar;

% HPBW

% Calcular el patrón de radiación
P = 2.*besselj(1,beta.*alpha.*sin(theta))./(beta.*alpha.*sin(theta));

% Encontrar el ángulo en el que la potencia es máxima
[max_intensity, max_index] = max(P);
max_angle = theta(max_index);

% Encontrar los ángulos en los cuales la potencia cae a la mitad del máximo
half_max_intensity = max_intensity / 2;
left_index = find(P(1:max_index) <= half_max_intensity, 1, 'last');
right_index = max_index + find(P(max_index:end) <= half_max_intensity, 1, 'first');

% Calcular el HPBW
hpbw = abs(theta(right_index) - theta(left_index));

disp(['HPBW:', num2str(hpbw)]);

% Graficar el patrón de radiación y las líneas para HPBW
figure;
plot(theta, P, 'b');
hold on;
plot([theta(left_index), theta(left_index)], [0, max_intensity], 'r--');
plot([theta(right_index), theta(right_index)], [0, max_intensity], 'r--');
hold off;
title('Patrón de Radiación 2 * J1(βαsin(θ)) / (βαsin(θ))');
xlabel('\theta');
ylabel('Intensidad');
grid on;
legend('Patrón de Radiación', 'Límites HPBW');

% Directividad

% Definir la función a integrar (ángulo sólido)
integrand = @(x, y) 2.*besselj(1,beta.*alpha.*sin(x))./(beta.*alpha.*sin(x)); 

% Definir los límites de integración para x (primer integral)
x_lower_limit = 0;
x_upper_limit = pi;

% Definir los límites de integración para y (segunda integral)
y_lower_limit = 0;
y_upper_limit = 2 * pi;

% Calcular la doble integral
result = integral2(integrand, x_lower_limit, x_upper_limit, y_lower_limit, y_upper_limit);

fprintf('Resultado de la doble integral: %.6f\n', result);

fprintf('Máxima Directividad: %.6f\n', 4*pi/result);
