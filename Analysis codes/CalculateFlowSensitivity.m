%CalculateFlowSensitivity.m
%Determine the sensitivity of flow to the spatial arrangement of landscape
%features, based on the Kaplan simulations. The final part of the code
%calculates the change in velocity that is attributable to each aspect of
%landscape configuration along a degradation trajectory. 
%Users should have first run CalculateOmegaHighFlow.m. 

load omegamodel.mat %Contains the variables ftot, Ctot, x1, x2, x3, x4, x5, omega_act, and omegacalc. This is created in the code OmegaModelingHighFlow.m.

P = 0.1:0.05:0.95; %Patch coverage. The last index is included for plotting only. Ignore values.
Hslough = 0.4:0.1:1.4; %Slough water depth, meters. The last index is included for plotting only. Ignore values.
DCI = 0.05:0.05:0.85; %Directional connectivity index
FD = 1.4:0.05:2; %Box-counting fractal dimension
Anis = 1:0.5:6; %anisotropy
S = 3e-5; %Water surface slope
n = 0.45; %Manning roughness, metric units

[p, h_slough, dci, anis] = ndgrid(P, Hslough, DCI, Anis); %Put these variables into a multidimensional matrix
fd = polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], p); %Regression fit for fd as a function of p. See SI for details.

%Reshape the landscape configuration matrices as linear arrays.
plin = reshape(p, numel(p),1);
h_sloughlin = reshape(h_slough, numel(h_slough),1);
dcilin = reshape(dci, numel(dci),1);
anislin = reshape(anis, numel(anis),1);
fdlin = reshape(fd, numel(fd),1);

%SET UP SINGLE-STANDARD-DEVIATION PERTURBATIONS
%Set up DCI perturbations
dci_high = min(1,dcilin+sqrt(sum((dcilin-polyval([10.741, -29.416, 24.237, 0.74982, -9.1403, 2.3807, 0.55588], plin)).^2)./numel(plin)));
dci_low = max(0,dcilin-sqrt(sum((dcilin-polyval([10.741, -29.416, 24.237, 0.74982, -9.1403, 2.3807, 0.55588], plin)).^2)./numel(plin)));

%Set up Hslough perturbations
hslough_high = h_sloughlin+std([0.4 0.5 0.6 0.7 0.8 0.9 1 1.2 1.4]);
hslough_low = h_sloughlin-std([0.4 0.5 0.6 0.7 0.8 0.9 1 1.2 1.4]);

%Set up fd perturbations
fd_high = min(2,fdlin+sqrt(sum((fdlin-polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], plin)).^2)./numel(plin)));
fd_low = fdlin-sqrt(sum((fdlin-polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], plin)).^2)./numel(plin));

%Set up p perturbations
p_high = min(1,plin+std([0.1 0.35 0.425 0.5 0.575 0.65 0.90]));
p_low = max(0,plin-std([0.1 0.35 0.425 0.5 0.575 0.65 0.90]));

%Set up anisotropy perturbations
anis_high = anislin+std([1 2 4 6]);
anis_low = max(1,anislin-std([1 2 4 6]));

%Solve for the baseline. 
h_ridge = max(h_slough-0.25,0);
omega = reshape(ftot(Ctot, [plin,dcilin,anislin,h_sloughlin,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_base = h_eff.^(2/3)*sqrt(S)/n;

%Now perturb DCI. 
omega = reshape(ftot(Ctot, [plin,dci_high,anislin,h_sloughlin,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
omega = reshape(ftot(Ctot, [plin,dci_low,anislin,h_sloughlin,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_dci = ((u_high-u_base)+(u_base-u_low))/2;
Du_dci(:,:,1,:) = u_high(:,:,1,:)-u_base(:,:,1,:); %Average of positive and negative perturbations.

%Now perturb h
h_ridge2 = reshape(hslough_high-0.25, size(h_slough));
omega = reshape(ftot(Ctot, [plin,dcilin,anislin,hslough_high,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge2.^omega + (1-p).*reshape(hslough_high, size(h_slough)).^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
h_ridge2 = reshape(max(0,hslough_high-0.25), size(h_slough));
omega = reshape(ftot(Ctot, [plin,dcilin,anislin,hslough_low,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge2.^omega + (1-p).*reshape(hslough_low, size(h_slough)).^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_h = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations.

%Now perturb fd
omega = reshape(ftot(Ctot, [plin,dcilin,anislin,h_sloughlin,fd_high]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
omega = reshape(ftot(Ctot, [plin,dcilin,anislin,h_sloughlin,fd_low]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_fd = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations.

%Now perturb p. 

fd2 = polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], p_high);
omega = reshape(ftot(Ctot, [p_high,dcilin,anislin,h_sloughlin,fd2]), size(h_slough)); 
h_eff = (reshape(p_high, size(p)).*h_ridge.^omega + (1-reshape(p_high, size(p))).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
fd2 = polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], p_low);
omega = reshape(ftot(Ctot, [p_low,dcilin,anislin,h_sloughlin,fd2]), size(h_slough)); 
h_eff = (reshape(p_low, size(p)).*h_ridge.^omega + (1-reshape(p_low, size(p))).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_p = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations.
Du_p(1:5,:,:,:)=u_high(1:5,:,:,:)-u_base(1:5,:,:,:); %Just use positive perturbations for these so as not to go out of range
Du_p(14:18,:,:,:) = u_base(14:18,:,:,:)-u_low(14:18,:,:,:); %Just use negative perturbations for these so as not to go out of range.

%Now perturb anisotropy. 
omega = reshape(ftot(Ctot, [plin,dcilin,anis_high,h_sloughlin,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
omega = reshape(ftot(Ctot, [plin,dcilin,anis_low,h_sloughlin,fdlin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_anis = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations.
Du_anis(:,:,:,1:3) = u_high(:,:,:,1:3)-u_base(:,:,:,1:3); %Just use positive perturbations for these so as not to go out of range.

%Plot results. Top row of plots is the well-connected landscape (see SI);
%bottom row is the poorly connected landscape.
figure
subplot(2,5,1)
pcolor(Hslough, P, Du_dci(:,:,15,8)), colorbar, title('DCI sensitivity')
subplot(2,5,2)
pcolor(Hslough, P, Du_h(:,:,15,8)), colorbar, title('h sensitivity')
subplot(2,5,3)
pcolor(Hslough, P, Du_fd(:,:,15,8)), colorbar, title('f_d sensitivity')
subplot(2,5,4)
pcolor(Hslough, P, real(Du_p(:,:,15,8))), colorbar, title('p sensitivity')
subplot(2,5,5)
pcolor(Hslough, P, Du_anis(:,:,15,8)), colorbar, title('e sensitivity')
subplot(2,5,6)
pcolor(Hslough, P, Du_dci(:,:,4,1)), colorbar, title('DCI sensitivity')
subplot(2,5,7)
pcolor(Hslough, P, Du_h(:,:,4,1)), colorbar, title('h sensitivity')
subplot(2,5,8)
pcolor(Hslough, P, Du_fd(:,:,4,1)), colorbar, title('f_d sensitivity')
subplot(2,5,9)
pcolor(Hslough, P, real(Du_p(:,:,4,1))), colorbar, title('p sensitivity')
subplot(2,5,10)
pcolor(Hslough, P, Du_anis(:,:,4,1)), colorbar, title('e sensitivity')


% NOW LOOK AT CHANGES IN VELOCITY AT TWO POINTS ALONG A DEGRADATION
% TRAJECTORY
%Base case landscape configuration parameters
dci = 0.6984;
h_ridge = 0.21-0.2;
h_slough = 0.72-0.2;
fd = 1.8544;
p = 0.4096;
anis = 3.5924;
omega = ftot(Ctot, [p,dci,anis,h_slough,fd]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_base = h_eff.^(2/3)*sqrt(S)/n;

%Perturbed case landscape configuration parameters
dci_p= 0.1523;
h_ridge_p = 0.20-0.2;
h_slough_p = 0.69-0.2;
fd_p = 1.7892;
p_p = 0.6025;
S_p = 2.5e-5;
anis_p = 2.8552;


labels = {'DCI', 'h', 'fd', 'p', 'S', 'e', 'DCI+e', 'p+fd', 'DCI+e+p+fd', 'h+S', 'All'};
u_pert = NaN(11,1);

%Perturbation due to DCI
omega = ftot(Ctot, [p,dci_p,anis,h_slough,fd]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_pert(1) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to h, water depth
omega = ftot(Ctot, [p,dci,anis,h_slough_p,fd]); 
h_eff = (p.*h_ridge_p.^omega + (1-p).*h_slough_p.^omega).^(1./omega);
u_pert(2) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to fractal dimension
omega = ftot(Ctot, [p,dci,anis,h_slough,fd_p]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_pert(3) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to patch coverage
omega = ftot(Ctot, [p_p,dci,anis,h_slough,fd]); 
h_eff = (p_p.*h_ridge.^omega + (1-p_p).*h_slough.^omega).^(1./omega);
u_pert(4) = h_eff.^(2/3)*sqrt(S)/n;

%Peturbation due to water surface slope
omega = ftot(Ctot, [p,dci,anis,h_slough,fd]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_pert(5) = h_eff.^(2/3)*sqrt(S_p)/n;

%Perturbation due to anisotropy
omega = ftot(Ctot, [p,dci,anis_p,h_slough,fd]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_pert(6) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to DCI + anisotropy
omega = ftot(Ctot, [p,dci_p,anis_p,h_slough,fd]); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_pert(7) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to patch cover + fractal dimension
omega = ftot(Ctot, [p_p,dci,anis,h_slough,fd_p]); 
h_eff = (p_p.*h_ridge.^omega + (1-p_p).*h_slough.^omega).^(1./omega);
u_pert(8) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to patch coverage + DCI + anisotropy + fractal dimension
omega = ftot(Ctot, [p_p,dci_p,anis_p,h_slough,fd_p]); 
h_eff = (p_p.*h_ridge.^omega + (1-p_p).*h_slough.^omega).^(1./omega);
u_pert(9) = h_eff.^(2/3)*sqrt(S)/n;

%Perturbation due to water depth+ water-surface slope
omega = ftot(Ctot, [p,dci,anis,h_slough_p,fd]); 
h_eff = (p.*h_ridge_p.^omega + (1-p).*h_slough_p.^omega).^(1./omega);
u_pert(10) = h_eff.^(2/3)*sqrt(S_p)/n;

%Perturbation due to all parameters
omega = ftot(Ctot, [p_p,dci_p,anis_p,h_slough_p,fd_p]); 
h_eff = (p_p.*h_ridge_p.^omega + (1-p_p).*h_slough_p.^omega).^(1./omega);
u_pert(11) = h_eff.^(2/3)*sqrt(S_p)/n;

%Plot the results as a bar chart.
to_plot = (u_pert([1:2, 5:10])-u_base)./(u_pert(11)-u_base);
plot_labels = {labels{1:2}, labels{5:10}};
figure, bar(to_plot)
set(gca, 'XTickLabels', plot_labels)
