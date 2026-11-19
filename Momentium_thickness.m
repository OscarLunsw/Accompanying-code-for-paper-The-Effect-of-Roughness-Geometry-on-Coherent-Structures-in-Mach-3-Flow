
%% Momentum Thickness Calculation

% Input data (column vectors)
% y   = wall-normal coordinate [m]
% rho = density profile [kg/m^3]
% u   = streamwise velocity profile [m/s]

% Freestream conditions
rho_e = rho_mean(end);    % Freestream density
U_e   = u(end);      % Freestream velocity

% Ensure column vectors
y   = y(:);
rho_mean = rho_mean(:);
u   = u(:);

% Density-weighted momentum thickness
integrand = (rho_mean .* u ./ (rho_e * U_e)) .* ...
            (1 - u ./ U_e);

% Numerical integration
theta = trapz(y, integrand);

% Display result
fprintf('Momentum thickness = %.6f m\n', theta);
fprintf('Momentum thickness = %.3f mm\n', theta*1000);
