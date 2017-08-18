% GenFig3B.m 
% This routine performs a comparison between a field-calibrated
% Kadlec model of flow through Everglades site U3 (by Choi and Harvey,
% Wetlands, 2014; hereafter abbreviated CH) an upscaled Manning model with
% optimized parameters, and a standard Manning model (through a homogeneous
% environment) with optimized roughness. It then generates the plot
% depicted as Fig. 3B in Larsen, Ma, and Kaplan, "How important is
% connectivity for surface-water fluxes? A generalized expression for flow
% through heterogeneous landscapes."

dchar = 0.08:0.01:0.52; %Characteristic water depth (h-tilde in eq. S20)
sigma = 0.56; %Sigma from eq. S20
S = 3e-5; %Water-surface slope
Psi = 3.3e7; %Psi from eq. S20
hs = 0.16:0.01:0.60; %Array of channel water depths


kadlec2 = Psi.*dchar.^sigma.*S./86400; %Eq. S20 with units conversion, yielding velocity in m/s

c_guess = [0.5 3e5/86400 -3 0.25 3e5/86400]; %Initial guess of parameter values for optimization. c(1) = patch coverage; c(2) = patch roughness; c(3) = slope of omega model; c(4) = intercept of omega model; c(5) = channel roughness 
model = @(c,x) ((c(1).*(c(2).*x).^(c(3).*(x-0.15)+c(4)) + (1-c(1)).*(c(5).*(x-0.15)).^(c(3).*(x-0.15)+c(4))).^(1./(c(3).*(x-0.15)+c(4)))).^(2/3).*(3e-5)^0.5; %Upscaled Manning model of velocity, in m/s
c_calc = lsqcurvefit(model, c_guess, hs, kadlec2); %Perform parameter optimization through nonlinear least-squares fitting
computed = model(c_calc, hs); %Use optimized parameter values to solve for computed velocity
figure, hold on, plot(hs, kadlec2, 'o') %Plot the CH model
plot(hs, computed, 'r-') %Plot the upscaled Manning model results
RMS1 = sqrt(sum((kadlec2-computed).^2)/length(kadlec2)); %RMS error
SSE = sum((kadlec2-computed).^2); %Sum of squared error of model
SST = sum((kadlec2-mean(kadlec2)).^2); %Sum of squares of CH model with respect to the mean
R2_1 = 1-SSE/SST; %R^2 value

manning = @(c,x) c.^(-1).*x.^(2/3).*sqrt(S); %Define the non-upscaled Manning model function
c_calc2 = lsqcurvefit(manning, 0.5, dchar, kadlec2); %Perform least-squares optimization to compute the non-upscaled n
RMS2 = sqrt(sum((kadlec2-manning(c_calc2, dchar)).^2)/length(kadlec2)); %RMS error
computed2 = manning(c_calc2, dchar); %Use optimized parameter values to solve for computed velocity
SSE2 = sum((kadlec2-computed2).^2); %Sum of squared error of model
R2_2 = 1-SSE2/SST; %R^2 value

%Generate plot.
hold on, plot(hs, computed2, 'k--')
legend('Field-calibrated Kadlec', 'Upscaled generalized force balance', 'Best-fit Manning')
xlabel('Channel depth, m')
ylabel('Velocity, m/s')
set(gca, 'FontSize', 14)