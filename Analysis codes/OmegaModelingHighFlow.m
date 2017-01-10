%OmegaModelingLowFlow.m
%This code is a companion to Larsen, Ma, and Kaplan, "How important is
%connectivity for surface-water fluxes? A generalized expression for flow
%through heterogeneous landscapes." It reads in the SWIFT2D simulation
%results, calculates the exact value of omega for each simulation, and then
%guides the user through finding an appropriate functional form for omega.
%This is the high-flow version of the code, for water depths greater than 15 cm
%above the top of the roughness elements. We recommend running this code
%section by section rather than all at once, as it is fairly interactive. 


%Load input files. See README for metadata
A = xlsread('Synthetic Landscapes input only.xlsx'); %Load water depths and Heff (effective water depth), m
Q = xlsread('Synthetic Landscapes Q.xlx'); %Load the discharge, m^3/s
Q_stats2 = reshape(Q, 14,30,4,7); %Dimensions are water depths, landscape realizations, anisotropy, patch coverage

%Read in the landscape configuration metrics
DCI = xlsread('Landscape metrics.xls');
DCIw_dir = DCI(:,3); %weighted, directed DCI
DCIw_undir = DCI(:,4); %weighted, undirected DCI
DCIu_dir = DCI(:,5); %unweighted, directed DCI
DCIu_undir = DCI(:,6); %unweighted, undirected DCI
fd = DCI(:,7); %fractal dimension
tort_dir = DCI(:,8); %directed spanning path tortuosity
tort_undir = DCI(:,9); %undirected spanning path tortuosity

Hslough = A(:,1); %same as depths, but repeated so that array has same length as number of domains
Heff=A(:,2:end); %Effective water depth, m
pct_ridge =[0.1 0.35 0.425 0.5 0.575 0.65 0.90]; %Patch coverage

% Inititalize variables for calculating omega
omega = NaN(size(A,1), length(pct_ridge));
c = NaN(1,3);
f = @(x,c) (c(1).*(max(c(2)-0.25,0))^x+(1-c(1)).*c(2).^x).^(1./x)-c(3); %The equation that will be solved to compute omega

% Solve for omega exactly
for ii = 1:size(A,1)
    for jj = 1:length(pct_ridge)
        c(1) = pct_ridge(jj);
        c(2) = Hslough(ii);
        c(3) = Heff(ii,jj);
        omega(ii,jj) = fzero(@(x) f(x,c), 1);
    end
end

omega_plot = reshape(omega, 14,30,4,length(pct_ridge)); %depths, reps, anisotropy, %ridge
depths = [0 0.1 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.2 1.4]; %slough depths
figure, % Shows plot of omega (y axis) vs. water depth for 4 different anisotropies and 7 different %ridge coverages

%Generate plot matrix of omega, where the four rows correspond to
%anisotropy classes (in descending order), and the seven columns correspond
%to the patch coverage. 
for ii = 1:4
    for jj = 1:7
        subplot(4,7, (ii-1)*7+jj)
        plot(depths(2:14), omega_plot(2:14,:,ii,jj))
        hold on, plot([0.4 0.4], [max(max(omega_plot(2:14,:,ii,jj))), min(min(omega_plot(2:14,:,ii,jj)))], 'k--') %dashed line at h = 0.4 m
        plot([0.25 0.25], [max(max(omega_plot(2:14,:,ii,jj))), min(min(omega_plot(2:14,:,ii,jj)))], 'k--') %dashed line at h = 0.25 m.
        axis tight
        set(gca, 'XLim', [0.1 1.4])
    end
end

%Reshape matrices into a more convenient form
Q_plot = reshape(Q, 14, 30, 4, length(pct_ridge));
omega_plot(Q_plot==0)=NaN; %Reset no-flow landscapes to NaNs.
Heff_plot = reshape(Heff, 14, 30, 4, length(pct_ridge));

%EVALUATE THE FIT OF DIFFERENT POLYNOMIALS
P = NaN(2,30,4,7);
mcalc = [];
mact = [];
for ii = 1:4 %This set of loops fits linear regressions to the high-water-depth part of the curves above
    for jj = 1:7
        for kk = 1:30
            x = transpose(depths(6:14))-0.25; 
            y = omega_plot(6:14,kk,ii,jj);
            P1 = polyfit(x(~isnan(y)),y(~isnan(y)),1);
            P(:,kk,ii,jj) = [P1'];%; P2'];
            mcalc = [mcalc; polyval(P1, x(~isnan(y)))];%Values calculated with the linear regression just computed above 
            mact = [mact; y(~isnan(y))]; %Actual values of omega
        end
    end
end
rsquared_tot = 1-sum((mcalc-mact).^2)/sum((mact-mean(mact)).^2); %R squared for the linear fits

%Now look at what the parameters vary with.
P_stats = shiftdim(P,1); %30 realizations,4 anisotropies,7 patch coverages,
P_stats = reshape(P_stats, 120,7,2); %replicates (30 realizations x 4 anisotropies), coverage, parameter

%Reshape landscape configuration matrices to a more convenient form
DCIu_dir_stats2 = reshape(DCIu_dir, 120, 7); 
DCIu_undir_stats2 = reshape(DCIu_undir, 120, 7);
DCIw_dir_stats2 = reshape(DCIw_dir, 120, 7); 
DCIw_undir_stats2 = reshape(DCIw_undir, 120, 7);
fd_stats2 = reshape(fd, 120, 7);
tort1stats2 = reshape(tort_dir, 120, 7);
tort2stats2 = reshape(tort_undir,120,7);
for ii = 1:2
    figure %Shows slope and intercept of the linear fits above as functions of fd, DCI undirected, and DCI directed
    for jj = 1:7
        subplot(7,7, jj)
        plot(DCIu_dir_stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'DCIu directed, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7,jj+7)
        plot(DCIu_undir_stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'DCIu undirected, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7, jj+14)
        plot(DCIw_dir_stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'DCIw directed, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7,jj+21)
        plot(DCIw_undir_stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'DCIw undirected, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7,jj+28)
        plot(fd_stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'fd, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7,jj+35)
        plot(tort1stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'tort dir, ', pct_ridge(jj), ' ridge coverage'))
        subplot(7,7,jj+42)
        plot(tort2stats2(:,jj), P_stats(:,jj,ii), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
        title(sprintf('%s%f%s', 'tort undir, ', pct_ridge(jj), ' ridge coverage'))
    end
end
    
%Now look at the parameters the above linear parameters vary with
%1. Slope
%THIS IS THE FINAL SLOPE MODEL!!
%x1 = patch coverage, x2 = anisotropy,x3 = fractal dimension
mact = reshape(P_stats(:,:,1), 840,1); %Actual value of slope
c_est = [0 0 -85.849 0 52.295 0 0 0 -0.881 -0.786 1 1]; %Vector of initial parameter estimates, based on Excel fitting
ftot = @(c,x) (c(1)+c(2).*x(:,3)+c(3).*x(:,1).^c(9)+c(4).*x(:,2)+c(5).*(x(:,3)).*x(:,1).^c(10)+c(6).*(x(:,3)).*(x(:,2))+c(7).*x(:,1).^c(11).*(x(:,2))+c(8).*x(:,1).^c(12).*(x(:,3)).*(x(:,2)));%Functional form of model
Ctot = lsqcurvefit(ftot, c_est, [p_mat anis_mat fd], mact); %Optimize parameter values
[beta,R,J,covB] = nlinfit([p_mat anis_mat fd], mact, ftot, Ctot); %Perform optimization usinga different routine, to get error values
SE = sqrt(diag(covB)); %Standard error of the estimates
mcalc = ftot(Ctot, [p_mat anis_mat fd]); %Estimated values of slope
t_est = beta'./SE; %t statistic
p_est = 1-tcdf(abs(t_est), length(mact)-1); %p-values
figure, plot(mcalc, mact, 'ko') %Predicted vs. actual values of slope.
rsquared_tot = 1-sum((mcalc-mact).^2)/sum((mact-mean(mact)).^2);%R2 for fit.

%2. Intercept
%2. Alternate way to solve for the intercept: Group by anisotropy ratio
%instead of by percent sawgrass:
DCI_dir_stats3 = reshape(shiftdim(reshape(DCIu_dir_stats2, 30, 4, 7), 2), 210, 4);
P_stats3 = reshape(shiftdim(reshape(P_stats(:,:,2), 30, 4, 7), 2), 210, 4);
p3 = NaN(4,2);
error_intercept = NaN(210,4);
figure %Plot the slope as a function of ln(DCI) and produce linear regressions of the relationship. The first plot corresponds to a DCI of 6; the last one to a DCI of one.
for ii = 1:4
    subplot(2,2,ii)
    plot(log(DCI_dir_stats3(:,ii)), P_stats3(:,ii), 'ko')
    p3(ii,:) = polyfit(log(DCI_dir_stats3(:,ii)), P_stats3(:,ii), 1);
    hold on
    x = -4.5:0.5:0;
    plot(x, polyval(p3(ii,:), x), 'k-')
    error_intercept(:,ii) = P_stats3(:,ii)-polyval(p3(ii,:), log(DCI_dir_stats3(:,ii)));
end

%Now look at how the error in the above relationship is distributed.
%First, break it out again by patch coverage, and by anisotropy ratio

%Reshape matrices for convenience:
error_intercept4 = reshape(error_intercept, 7, 30, 4);
DCI_dir_stats4 = reshape(DCI_dir_stats3, 7, 30, 4);
fd_stats4 = reshape(reshape(shiftdim(reshape(fd_stats2, 30, 4, 7), 2), 210, 4), 7, 30, 4);
for ii = 1:7 %Loop over patch coverages
    for jj = 1:4 %Loop over anisotropies
        figure(50), subplot(4,7, sub2ind([7 4], ii, jj))
        plot(DCI_dir_stats4(ii,:,jj), error_intercept4(ii,:,jj), 'b*')
        axis tight
        figure(51), subplot(4,7, sub2ind([7 4], ii, jj))
        plot(fd_stats4(ii,:,jj), error_intercept4(ii,:,jj), 'k*')
        axis tight
    end
end

%Nonlinear possibility for the intercept. However,the final intercept model
%that we used was fitted as a multiple linear regression on JMP statistical
%software. The whole-omega model below incorporates it accurately.
%x1 = patch coverage, x2 = DCI, x3 = anisotropy
ftot = @(c,x) (c(1).*x(:,3).^c(2)+c(3)).*log(x(:,2))+c(4).*log(x(:,3))+c(5)+c(6).*x(:,1); %Functional form of this model
c_est = [0.8586 -0.511 -1.18 0.7651 -0.6816+0.6241 -2.3919]; %Vector of initial parameter estimates
anis_mat = reshape(repmat([6*ones(30,1); 4*ones(30,1); 2*ones(30,1); ones(30,1)], 1, 7), 840, 1);
p_mat = reshape(repmat(pct_ridge, 120, 1), 840, 1);
Ctot = lsqcurvefit(ftot, c_est, [p_mat, DCIw_undir, anis_mat], reshape(P_stats(:,:,2), 840,1)); %Optimize model parameters
[beta, R,J,covB] = nlinfit([p_mat, DCIw_undir, anis_mat], reshape(P_stats(:,:,2), 840,1), ftot,Ctot); %Do it another way to get error values
sqrt(diag(covB)) %the standard errors
bcalc = ftot(Ctot, [p_mat, DCIw_undir, anis_mat]); %Calculated intercept
bact = reshape(P_stats(:,:,2), 840,1); %Actual intercept
figure, plot(bcalc, bact, 'ko')
rsquared_tot = 1-sum((bcalc-bact).^2)/sum((bact-mean(bact)).^2);

%THIS IS THE FINAL COMBINED MODEL
%x1 = patch coverage, x2 = DCI, x3 = anisotropy, x4 = depth, x5 = fd
ftot = @(c,x) c(1)+c(2).*log(x(:,3))+c(3).*log(x(:,2))+c(4).*x(:,1)+c(5).*x(:,5)+c(6).*log(x(:,3)).*log(x(:,2))+c(7).*(x(:,1)-0.5).*(x(:,5)-1.77072) + (c(8)+c(9).*x(:,5)+c(10).*x(:,1).^c(16)+c(11).*x(:,3)+c(12).*(x(:,5)).*x(:,1).^c(17)+c(13).*(x(:,5)).*(x(:,3))+c(14).*x(:,1).^c(18).*(x(:,3))+c(15).*x(:,1).^c(19).*(x(:,5)).*(x(:,3))).*(x(:,4)-0.25); %functional form
c_est = [-9.829519, 0.5319014, 0.1922771, 2.3445588, 4.7571396, -0.092516, -4.061362...
     523.466584042198,-283.857604078859,-543.607919148115,0.350606730526567,299.714896502203,-0.160747417173356,-39.1990662152862,32.3955461844169,-0.294855296430155,-0.274958352840607,14.6018788489946,16.5897640670315] ; %Vector of parameter estimates
     
%Rearrange landscape configuration parameters into a more convenient form,
%in accordance with x1-5 nomenclature above.
x1 = shiftdim(repmat(pct_ridge', 1, 14, 30, 4), 1);
x2 = repmat(shiftdim(reshape(DCIu_dir_stats2, 30, 4, 7), -1),14, 1, 1, 1); 
x3 = shiftdim(repmat([6;4;2;1], 1, 7, 14, 30),2);
x4 = repmat(depths', 1, 30, 4, 7);
x5 = repmat(shiftdim(reshape(fd_stats2, 30, 4, 7), -1),14, 1, 1, 1); 

%Now these are in the dimensions of omega_plot: depths x reps x anisotropy x %ridge

x1 = reshape(x1(6:14, :, :, :), 9*30*4*7,1);
x2 = reshape(x2(6:14, :, :, :), 9*30*4*7,1);
x3 = reshape(x3(6:14, :, :, :), 9*30*4*7,1);
x4 = reshape(x4(6:14, :, :, :), 9*30*4*7,1);
x5 = reshape(x5(6:14, :, :, :), 9*30*4*7,1);

% % %See how we've done. Pretty good.
omega_check = ftot(c_est, [x1 x2 x3 x4 x5]);
figure, plot(omega_check, reshape(omega_plot(6:14, :, :, :), 9*30*4*7,1), 'k*')
rsquared_tot = 1-sum((omega_check-omega_act).^2)/sum((omega_act-mean(omega_act)).^2);


omega_act = reshape(omega_plot(6:14, :, :, :), 9*30*4*7,1);
Ctot = lsqcurvefit(ftot, c_est, [x1 x2 x3 x4 x5], omega_act);
omegacalc = ftot(Ctot, [x1 x2 x3 x4 x5]);
figure, plot(omegacalc, omega_act, 'ko')
rsquared_tot = 1-sum((omegacalc-omega_act).^2)/sum((omega_act-mean(omega_act)).^2);
save('omegamodel.mat', 'ftot', 'Ctot', 'x1', 'x2', 'x3', 'x4', 'x5', 'omega_act', 'omegacalc')


%Now check the fit of Q...
% Predicted Q
W = 2000;
n = 0.45;
S = 3e-5;
heff = ((1-x1).*x4.^omegacalc+x1.*(max(0,x4-.25)).^omegacalc).^(1./omegacalc);
heff_act = ((1-x1).*x4.^omega_act+x1.*(max(0,x4-.25)).^omega_act).^(1./omega_act);

Q_stats = reshape(Q_stats2(6:14, :, :, :), 9*30*4*7,1);
Q_check = W./n.*heff.^(5/3).*sqrt(S);
u_stats = Q_stats/W./heff_act;
u_check = Q_check/W./heff;
figure, plot(u_check,u_stats, 'k*')
figure, plot(u_stats, u_check-u_stats, 'k*')

