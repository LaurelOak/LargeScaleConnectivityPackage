% OmegaTheoretical.m
% This routine generates computations of omega for theoretical solutions of
% flow over a landscape with patches in parallel and a landscape with
% patches in series. This code was used to generate Figs. 1A-B in Larsen,
% Ma, and Kaplan, "How important is connectivity for surface-water fluxes?
% A generalized expression for flow through heterogeneous landscapes." To
% generate those figures, values of the lumped roughness coefficients and
% patch coverage were selected manually.

%% DEFINE CONSTANTS AND SOLVER FUNCTION
n = 0.45; %Manning roughness
H = 0.3:0.05:1.4; %Array of water depths in channel
g = 9.8; %Gravitational constant, m/s^2
delta_z = 0.25; %Patch relief, m
F = @(x,c) ((1-c(1))*c(2)^x+c(1)*(c(2)-0.25)^x)^(1/x)-c(3); %Solver function for omega. c(1) is patch coverage; c(2) is channel water depth; c(3) is Heff.


%% %% PATCHES IN PARALLEL: Fig. 1A
Heff = zeros(size(H)); %Initialize matrix of H_eff
omega = Heff; %Initialize omega
Wc = 400; %Width of channels, m
beta = 3.96; %Lumped drag coefficient
p = 0.1; %Patch coverage
for ii = 1:length(H)
    Heff(ii) = (n*(p*(H(ii)-delta_z).*sqrt(2*g*(H(ii)-delta_z).^(4/3)/beta)+(1-p)*H(ii).*sqrt(2*g*H(ii).^(4/3)/(beta+2*beta*delta_z/Wc)))).^(3/5); %Eq. S9
    c(1) = p; %Input to function to solve for omega
    c(2) = H(ii); %Input to function to solve for omega
    c(3) = Heff(ii); %Input to function to solve for omega
    omega(ii) = fzero(@(x) F(x,c), -1); %Solve for omega
end

% figure
hold on
plot(H, omega, 'k-', 'LineWidth', 0.5, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k') %Plot omega over a range of water depths. To generate Fig. 1A, I would change beta and p manually and update the plot.


%% PATCHES IN SERIES: Fig. 1B
Heff = zeros(size(H)); %Initialize matrix of H_eff
omega = Heff; %Initialize omega
L = 400; %Length of ridges, m
beta2 = 0.205; %Lumped bed drag coefficient
beta1 = 2; %Lumped patch form drag coefficient
Ld = 4000; %Length of domain, m
p = 0.1; %Patch coverage

for ii = 1:length(H) %Solve over a range of water depths
    Heff(ii) = (n*(H(ii)-delta_z)*sqrt((p*(H(ii)-delta_z)+(1-p)*H(ii))/(beta1*delta_z/L^2/H(ii)^(1/3)*p*Ld*(p+(1-p)*(H(ii)-delta_z)/H(ii))^2+beta2/(H(ii)-delta_z)^(1/3)*p+(1-p)*beta2*(H(ii)-delta_z)^2/H(ii)^(7/3))))^(3/5); %Equation S12
    c(1) = p; %Input to function to solve for omega
    c(2) = H(ii); %Input to function to solve for omega
    c(3) = Heff(ii); %Input to function to solve for omega
    omega(ii) = fzero(@(x) F(x,c), -1); %Solve for omega
end

% figure
hold on
plot(H, omega, 'r-', 'LineWidth', 0.5) %Plot omega over a range of water depths. To generate Fig. 1B, I would change beta and p manually and update the plot.
