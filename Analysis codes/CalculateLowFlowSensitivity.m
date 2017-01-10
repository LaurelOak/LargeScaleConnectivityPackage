%CalculateLowFlowSensitivity.m Determine the sensitivity of flow to the
%spatial arrangement of landscape features, based on the Kaplan
%simulations. This version of the code applies to water depths equal to or
%lower than the top of the roughness elements.

%CHOOSE A MODEL--h = zp, or h < zp. The other one should be commented.
%h = zp
ftot = @(c,x) c(1) + c(2).*x(:,2) + c(3).*x(:,1).^c(4) + c(5).*(x(:,1)-0.35192).^c(6).*(x(:,2)-0.47181)+c(7).*log(x(:,3))+c(8).*log(x(:,3)).*(x(:,1)-0.35192).^c(9); %Final version of the model
Ctot = [0.0585215, 0.2919471, 0.4051182, 1, -1.164312, 1, 0.1400144, -0.211838, 1]; %Final vector of coefficients

% %h<zp
% ftot = @(c,x) c(1) + c(2).*x(:,1) + c(3).*x(:,2) + c(4).*(x(:,1)-0.35192).*(x(:,2)-0.47181)+c(5).*log(x(:,3))+c(6).*log(x(:,3)).*(x(:,1)-0.35192); %Final version of the model
% Ctot = [-0.080711, 0.5475127, 0.4546676, -1.521892, 0.0840915, -0.369997]; %Final vector of coefficients


P = 0.1:0.05:0.95; %Patch coverage. The last index is included for plotting only. Ignore values.
Hslough = 0.20; %Channel water depth, meters. The last index is included for plotting only. Ignore values.
DCI = 0.05:0.05:0.85; %Directional connectivity index
FD = 1.4:0.05:2; %Box-counting fractal dimension
Anis = 1:0.5:6; %Anisotropy
S = 3e-5; %Water-surface slope
n = 0.45; %Manning roughness value

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
omega = reshape(ftot(Ctot, [plin,dcilin,anislin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_base = h_eff.^(2/3)*sqrt(S)/n;

%Now perturb DCI. 
omega = reshape(ftot(Ctot, [plin,dci_high,anislin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
omega = reshape(ftot(Ctot, [plin,dci_low,anislin]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_dci = ((u_high-u_base)+(u_base-u_low))/2;%Average of positive and negative perturbations.
Du_dci(:,:,1,:) = u_high(:,:,1,:)-u_base(:,:,1,:); %Just use positive perturbations for these so as not to go out of range

%Now perturb h
h_ridge2 = reshape(hslough_high-0.25, size(h_slough));
omega = reshape(ftot(Ctot, [plin,dcilin,anislin]), size(h_slough)); 
h_eff = (p.*h_ridge2.^omega + (1-p).*reshape(hslough_high, size(h_slough)).^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
h_ridge2 = reshape(max(0,hslough_low-0.25), size(h_slough));
omega = reshape(ftot(Ctot, [plin,dcilin,anislin]), size(h_slough)); 
h_eff = (p.*h_ridge2.^omega + (1-p).*reshape(hslough_low, size(h_slough)).^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_h = ((u_high-u_base)); %AT these low water depths, only use the positive perturbation.

%Now perturb p. 

fd2 = polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], p_high);
omega = reshape(ftot(Ctot, [p_high,dcilin,anislin]), size(h_slough)); 
h_eff = (reshape(p_high, size(p)).*h_ridge.^omega + (1-reshape(p_high, size(p))).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
fd2 = polyval([-31.442, 91.055, -104.79, 60.124, -17.881, 2.1873, 1.8562], p_low);
omega = reshape(ftot(Ctot, [p_low,dcilin,anislin]), size(h_slough)); 
h_eff = (reshape(p_low, size(p)).*h_ridge.^omega + (1-reshape(p_low, size(p))).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_p = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations.
Du_p(1:5,:,:,:)=u_high(1:5,:,:,:)-u_base(1:5,:,:,:); %Just use positive perturbations for these so as not to go out of range
Du_p(14:18,:,:,:) = u_base(14:18,:,:,:)-u_low(14:18,:,:,:);%Just use negative perturbations for these so as not to go out of range

%Now perturb anisotropy. 
omega = reshape(ftot(Ctot, [plin,dcilin,anis_high]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_high = h_eff.^(2/3)*sqrt(S)/n;
omega = reshape(ftot(Ctot, [plin,dcilin,anis_low]), size(h_slough)); 
h_eff = (p.*h_ridge.^omega + (1-p).*h_slough.^omega).^(1./omega);
u_low = h_eff.^(2/3)*sqrt(S)/n;
Du_anis = ((u_high-u_base)+(u_base-u_low))/2; %Average of positive and negative perturbations
Du_anis(:,:,:,1:3) = u_high(:,:,:,1:3)-u_base(:,:,:,1:3); %Just use positive perturbations for these so as not to go out of range

%Plot results. Top row of plots is the well-connected landscape (see SI);
%bottom row is the poorly connected landscape.
figure
subplot(2,5,1)
pcolor(Anis, P, squeeze(Du_dci(:,:,15,:))), colorbar, title('DCI sensitivity')
subplot(2,5,2)
pcolor(Anis, P, squeeze(Du_h(:,:,15,:))), colorbar, title('h sensitivity')
subplot(2,5,4)
pcolor(Anis, P, squeeze(real(Du_p(:,:,15,:)))), colorbar, title('p sensitivity')
subplot(2,5,5)
pcolor(Anis, P, squeeze(Du_anis(:,:,15,:))), colorbar, title('e sensitivity')
subplot(2,5,6)
pcolor(Anis, P, squeeze(Du_dci(:,:,4,:))), colorbar, title('DCI sensitivity')
subplot(2,5,7)
pcolor(Anis, P, squeeze(Du_h(:,:,4,:))), colorbar, title('h sensitivity')
subplot(2,5,9)
pcolor(Anis, P, squeeze(real(Du_p(:,:,4,:)))), colorbar, title('p sensitivity')
subplot(2,5,10)
pcolor(Anis, P, squeeze(Du_anis(:,:,4,:))), colorbar, title('e sensitivity')
