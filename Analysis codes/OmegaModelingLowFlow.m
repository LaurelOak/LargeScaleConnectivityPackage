%OmegaModelingLowFlow.m
%This code is a companion to Larsen, Ma, and Kaplan, "How important is
%connectivity for surface-water fluxes? A generalized expression for flow
%through heterogeneous landscapes." It reads in the SWIFT2D simulation
%results, calculates the exact value of omega for each simulation, and then
%guides the user through finding an appropriate functional form for omega.
%This is the low-flow version of the code, for water depths less than 15 cm
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
        plot([0.25 0.25], [max(max(omega_plot(2:14,:,ii,jj))), min(min(omega_plot(2:14,:,ii,jj)))], 'k--') %dashed line at h = 0.25 m
        axis tight
        set(gca, 'XLim', [0.1 1.4])
    end
end
Q_plot = reshape(Q, 14, 30, 4, length(pct_ridge));
omega_plot(Q_plot==0)=NaN; %Reset omega values for no-flow landscapes to NaN
Heff_plot = reshape(Heff, 14, 30, 4, length(pct_ridge));

this_omega = squeeze(omega_plot(4,:,:,:)); %Water level at top of ridges

%Put landscape configuration matrices into a more convenient form
DCIu_dir_stats2 = reshape(DCIu_dir, 120, 7); 
DCIu_undir_stats2 = reshape(DCIu_undir, 120, 7);
DCIw_dir_stats2 = reshape(DCIw_dir, 120, 7); 
DCIw_undir_stats2 = reshape(DCIw_undir, 120, 7);
fd_stats2 = reshape(fd, 120, 7);
tort1stats2 = reshape(tort_dir, 120, 7);
tort2stats2 = reshape(tort_undir,120,7);
    
%Plot omega as a function of DCI, grouped by pct coverage, and fit a linear
%regression through it.
P_stats = NaN(2,7);
resids = NaN(120,7);
figure
for ii=1:7
    x = DCIu_dir_stats2(:,ii);
    y = reshape(this_omega, 120, 7);
    y = y(:,ii);
    x = x(~isnan(y));
    y = y(~isnan(y));
    P_stats(:,ii) = polyfit(x,y,1);
    subplot(1,7,ii), plot(x,y,'o', min(x):0.01:max(x), polyval(P_stats(:,ii), min(x):0.01:max(x)), 'k-')
    x = DCIu_dir_stats2(:,ii);
    y = reshape(this_omega, 120, 7);
    y = y(:,ii);
    resids(:,ii) = y-polyval(P_stats(:,ii), x); %Residuals of the regression
end


%Plot omega as function of anisotropy, grouped by patch coverage
P_stats = NaN(2,7);
resids = NaN(120,7);
figure
for ii=1:7
    x = reshape(anis_mat, 120,7); 
    x = (x(:,ii));
    y = reshape(this_omega, 120, 7);
    y = y(:,ii);
    x = log(x(~isnan(y)));
    y = y(~isnan(y));
    P_stats(:,ii) = polyfit(x,y,1);
    subplot(1,7,ii), plot(x,y,'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.0), hold on, plot( min(x):0.01:max(x), polyval(P_stats(:,ii), min(x):0.01:max(x)), 'k-')
    set(gca, 'FontSize', 14), xlabel('ln(e)')
end

%Fit parameters to one possible model with just DCI and patch coverage.
%Ultimately this one wasn't used. A stepwise multiple linear regression,
%performed with JMP statistical software, was the final model used. The
%model and parameter estimates for the model can be found in the first few
%lines of the accompanying sensitivity analysis code.
%x(1) = patch coverage, x(2) = DCI 
ftot = @(c,x) c(1).*x(:,1).^c(2).*x(:,2)+c(3).*x(:,1).^c(4)+c(5);
c_est = [0.3854 -1.039 -0.0955 -1.546 0.05];
y = reshape(this_omega, 840,1);
pct_ridge1 = reshape(repmat(pct_ridge, 120, 1), 840,1);
x = [pct_ridge1, DCIu_dir];
x = x(~isnan(y),:);
y = y(~isnan(y));
Ctot = lsqcurvefit(ftot,c_est,x,y); %Optimize parameters
figure, plot(ftot(Ctot,x),y,'ko') %Plot predicted vs actual
these_resids = reshape(this_omega,840,1)-ftot(Ctot, [pct_ridge1, DCIu_dir]); %Residuals
calced = ftot(Ctot,x); %Predicted values of omega
rsquared_tot = 1-sum((calced-y).^2)/sum((y-mean(y)).^2);
[beta,R,J,covB] = nlinfit(x, y, ftot, Ctot); %Solve using nlinfit to get error estimates
SE = sqrt(diag(covB)); %Standard error of the estimates
t_est = beta'./SE; %T-statistics
p_est = 1-tcdf(abs(t_est), length(mact)-1); %P-values



%2. LOW LOW WATER FITTING (h<zp)
A1 = reshape(squeeze(omega_plot(3,:,:,:)), 840, 1); %Water level at h = 0.2
A2 = reshape(squeeze(omega_plot(2,:,:,:)), 840, 1); %Water level at h = 0.1
LFO = NaN(size(A2)); %Initialize the low-flow omega
%CAlculate an average low-flow omega to model
for ii = 1:840
    if isnan(A2(ii))
        LFO(ii) = A1(ii);
    else LFO(ii) = (A1(ii)+A2(ii))/2;
    end
end
LFO(LFO>1.6)=NaN; %There was one outlier with an anomalously high value. Remove it.

%Below, plot the low-flow omega as a function of the DCI, and fit a linear
%regression through the relationship.
lfo = reshape(LFO, 120, 7);
figure
P_stats = NaN(2,7);
for ii = 1:7
    subplot(1,7,ii)
    x = DCIw_dir_stats2(:,ii);
    y = (lfo(:,ii));
    x = x(~isnan(y)); %Remove no-flow NaNs.
    y = y(~isnan(y));
    plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
    xlabel('DCI')
    title(sprintf('%s%d', 'p = ', pct_ridge(ii)))
    P_stats(:,ii) = polyfit(x, y, 1);
    hold on
    plot([min(x), max(x)], polyval(P_stats(:,ii), [min(x), max(x)]), 'k-')
end

%Plot the low-flow omega as a function of anisotropy
figure
for ii=1:7
    x = reshape(anis_mat, 120,7); 
    x = (x(:,ii));
    y = lfo(:,ii);
    x = log(x(~isnan(y)));
    y = y(~isnan(y));
    P_stats(:,ii) = polyfit(x,y,1);
    subplot(1,7,ii), plot(x,y,'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.0), hold on, plot( min(x):0.01:max(x), polyval(P_stats(:,ii), min(x):0.01:max(x)), 'k-')
    set(gca, 'FontSize', 14), xlabel('ln(e)')
end

%Generate and fit one possible model that is only a function of the DCI.
%Ultimately, this one was not used. x1 = patch coverage, x2 = DCI
ftot = @(c,x) (c(1).*x(:,1)+c(2)).*exp(c(3).*x(:,1).^c(4).*x(:,2));
c_est = [0.6502 -0.0692 0.8418 -0.992];
y = LFO;
x = [pct_ridge1, DCIw_dir];
x = x(~isnan(y),:);
y = y(~isnan(y));
Ctot = lsqcurvefit(ftot,c_est,x,y);
[beta,R,J,covB] = nlinfit(x, y, ftot, Ctot);
SE = sqrt(diag(covB)); %Standard error of the estimates
t_est = beta'./SE;%t-statistic of estimates
p_est = 1-tcdf(abs(t_est), length(mact)-1); %p-values associated with estimates

figure, plot(ftot(Ctot,x),y,'ko') %plot predicted vs actual omega
calced = ftot(Ctot,x); %predicted version of omega
rsquared_tot = 1-sum((calced-y).^2)/sum((y-mean(y)).^2);
these_resids = LFO-ftot(Ctot, [pct_ridge1, DCIw_dir]); %residuals

%Another alternative, which combines the model above with a multiple linear
%regression on patch coverage and anisotropy. x1 = patch coverage, x2 =
%DCI, x3 = anisotropy. 
%Ultimately, the model that was used was a multiple linear regression fit
%on JMP statistical software. It can be found in the first few lines of the
%sensitivity analysis code.
ftot = @(c,x) (c(1).*x(:,1)+c(2)).*exp(c(3).*x(:,1).^c(4).*x(:,2)) + c(5).*x(:,1) + c(6).(x(:,3)) + c(7).*x(:,1).*(x(:,3));
c_est = [0.6502 -0.0692 0.8418 -0.992 0 0 0]; %Initial parameter estimates
y = LFO;
x = [pct_ridge1, DCIw_dir, anis_mat];
x = x(~isnan(y),:);
y = y(~isnan(y));
Ctot = lsqcurvefit(ftot,c_est,x,y); %Optimize parameters
[beta,R,J,covB] = nlinfit(x, y, ftot, Ctot); %Optimize parameters another way, in order to get errors.
SE = sqrt(diag(covB)); %Standard error of the estimates
t_est = beta'./SE; %t-statistic
p_est = 1-tcdf(abs(t_est), length(y)-1); %p-value

figure, plot(ftot(beta,x),y,'ko')%Predicted vs actual omega
calced = ftot(beta,x); %Predicted omega
rsquared_tot = 1-sum((calced-y).^2)/sum((y-mean(y)).^2);
these_resids = LFO-ftot(beta, [pct_ridge1, DCIw_dir]); %Residuals

