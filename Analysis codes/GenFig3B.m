% GenFig3B.m 
% This routine performs a comparison between a field-calibrated
% Kadlec model of flow through Everglades site U3 (by Choi and Harvey,
% Wetlands, 2014; hereafter abbreviated CH) an upscaled Manning model with
% optimized roughness parameters, and a standard Manning model (through a 
% homogeneous environment) with optimized roughness. It then generates the 
% plot depicted as Fig. 3B in Larsen, Ma, and Kaplan, "How important is
% connectivity for surface-water fluxes? A generalized expression for flow
% through heterogeneous landscapes."

dchar = 0.08:0.01:0.52; %Characteristic water depth (h-tilde in eq. S20) 
sigma = 0.56; %Sigma from eq. S20
S = 3e-5; %Water-surface slope
Psi = 3.3e7; %Psi from eq. S20
hs = 0.16:0.01:0.60; %Array of channel water depths 

DCI = 0.1016; %unweighted, directed directional connectivity index for landscape PSU 37
DCI_w = 0.0265; %weighted, directed directional connectivity index for landscape PSU 37
anis = 1.0246; %anisotropy for landscape PSU 37
pct_ridge = 0.856; %patch coverage for landscape PSU 37
fd = 1.3755; %box-counting fractal dimension for landscape PSU 37

%Based on the measured landscape parameters above, calculate the slope and
%intercept of the omega-vs-water-depth relation, in accordance with Table
%S1 and equation S12:

m = 523-284*fd-543*pct_ridge^-0.29+0.35*anis+300*fd*pct_ridge^-0.28-0.16*anis*fd-39*pct_ridge^14.6*anis+32*pct_ridge^16.6*anis*fd; %Computed slope
b = -9.83+0.53*log(anis)+2.34*pct_ridge+4.76*fd-4.06*(pct_ridge-0.5)*(fd-1.77)+.19*log(DCI)-0.09*log(DCI)*log(anis); %Computed intercept

omega = m*(hs-0.11)+b; %Array of omega over water depths. 0.11 = zp, measured vertical relief of patches, from Choi and Harvey 2014.
upper_lim_omega_interp = omega(hs==0.26); %Interpolate transitional values between this and omega_low
omega_low = 0.14+0.22*DCI_w+0.25*pct_ridge-0.65*(pct_ridge-0.4)*(DCI_w-0.38)+0.14*log(anis)-0.14*log(anis)*(pct_ridge-0.40); % Modeled omega for h = zp (Table S2).
omega(hs<0.26)= interp1([0.16, 0.26], [omega_low, upper_lim_omega_interp], hs(hs<0.26)); %Linearly interpolate omega for transitional water depths

kadlec2 = Psi.*dchar.^sigma.*S./86400; %Eq. S20 with units conversion, yielding velocity in m/s

model = @(c,x) ((pct_ridge.*(c(1).*x).^(m.*(x-0.11)+b) + (1-pct_ridge).*(c(2).*(x-0.11)).^(m.*(x-0.11)+b)).^(1./(m.*(x-0.11)+b))).^(2/3).*(3e-5)^0.5; %Upscaled Manning model of velocity, in m/s. c(1) = K_patch; c(2) = K_channel
c_est = [2.8614, 4.6525]; %Initial guess for K ([patch, channel])
c_calc = lsqcurvefit(model, c_est, hs, kadlec2); %Optimize the two K-values.
computed = model(c_calc, hs); %Generate estimate of velocity based on the optimized K-values.
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