%% HW 2
close all;
load hw2.mat;
fileName = fopen('hw2.txt','w');
%% Problem 1.a

% Separate data into x and y
xData = data1(:,1);
yData = data1(:,2:end);
yAvg = zeros(length(xData), 1);
yStd = zeros(length(xData), 1);
ySEM = zeros(length(xData), 1);

% Find the mean and std at each x value
for i = 1:length(xData)
    yAvg(i) = mean(yData(i,:));
    yStd(i) = std(yData(i,:));
    ySEM(i) = yStd(i) / sqrt(length(xData));
end
fprintf(fileName,'x value | Mean | Uncertainty | Std  \n');
fprintf(fileName,'%d \t %d \t %d \t %d \n',xData,yAvg,ySEM,yStd);

% Find chi^2 min
fit = @(p, xData) p(1) .* xData .^p(2); % Power fit
chi2 = @(p) sum( ((yAvg - fit(p,xData)) ./ ySEM).^2); % Chi^2 for polynomial

% Find parameters to minimize chi^2
pmin = fminsearch(chi2, [0 0]);
chi2_min = chi2(pmin);
chi2_reduced = chi2_min / (length(xData) - 2);
fprintf(fileName,'A = %f, p = %f \n', pmin(1), pmin(2));
fprintf(fileName,'Chi^2_min = %f, Chi^2_reduced = %f \n', chi2_min, chi2_reduced);

% Find uncertainty in A and p using Chi^2 + 1
A = pmin(1);
p = pmin(2);
rA = .5;
rP = rA;
ARange = linspace(A - rA, A + rA, 100);
pRange = linspace(p - rP, p + rP, 100);

% 1. Generate a grid to search parameter pairs for min Chi^2 + 1
[ASpace, pSpace] = ndgrid(ARange, pRange); 
chiSquareGrid = zeros(size(ASpace));

% 2. In the grid, compute Chi^2 for every parameter combination
for i = 1:size(ASpace,1)
    for j = 1:size(ASpace,2)
        chiSquareGrid(i,j) = chi2([ASpace(i,j) pSpace(i,j)]);
    end
end    

% 3. Find the indices where chi^2 <= chi2_min + 1
[indexA, indexP] = find(chiSquareGrid <= chi2_min + 1);

% 4. Find the min and max of these parameters
dAMin = min(ARange(indexA));
dAMax = max(ARange(indexA));
dpMin = min(pRange(indexP));
dpMax = max(pRange(indexP));

% 5. Compute the uncertainties!
dA = (dAMin - dAMax)/2;
dP = (dpMin - dpMax)/2;
fprintf(fileName,'A Uncertainty: %f \t B Uncertainty: %f \n', dA, dP);

% Plot to check!
% figure
% hold on
% contourf(ASpace, pSpace, chiSquareGrid, 50);
% scatter(ARange(indexA), pRange(indexP), 'mo', 'MarkerFaceColor', 'Magenta');
% plot(pmin(1), pmin(2), 'kx', 'MarkerSize', 15);
% xlabel('A'); ylabel('b'); ch = colorbar; ylabel(ch, sprintf('\\chi^2'));
% title(sprintf('\\chi^2 Parameter Space'));

%% Problem 1.b
% OLS
% 1. Linearize data
x_lin = log(xData);
y_lin = log(yAvg);

% 2. Get y_tilde = A_tilde + p*x_tilde and do OLS
x_ols = [ones(length(x_lin),1),x_lin];
p_ols = x_ols \ y_lin;

% 3. Rescale p_ols back to get A and p
AVal_ols = exp(p_ols(1));
pVal_ols = p_ols(2);
Rsq_ols = 1 - sum((yAvg - fit([AVal_ols,pVal_ols], xData)).^2)/...
    sum((yAvg - mean(yAvg)).^2);
fprintf(fileName,'OLS Fit: A = %f, p = %f, Rsq = %f\n', ...
    AVal_ols, pVal_ols, Rsq_ols);

%% Problem 1.c
% WLS
% 1. Weight x and y based on their uncertainty
% NOTE: The log comes from linearization
x_wls = x_ols ./ ySEM;
y_wls = y_lin ./ ySEM;

% 2. Go through the least squares method
p_wls = x_wls \ y_wls;
AVal_wls = exp(p_wls(1));
pVal_wls = p_wls(2);
Rsq_wls = 1 - sum((yAvg - fit([AVal_wls,pVal_wls], xData)).^2)/...
    sum((yAvg - mean(yAvg)).^2);

% 3. Find the uncertainties in A and p
cov = inv(x_wls'*x_wls);
sigma_param = sqrt(diag(cov));
dA_wls = sigma_param(1);
dp_wls = sigma_param(2);

fprintf(fileName,'WLS Fit: A = %f, p = %f, Rsq = %f\n', ...
    AVal_wls, pVal_wls, Rsq_wls); 
fprintf(fileName,'WLS Uncertainties: dA = %f, dp = %f\n', ...
    dA_wls, dp_wls);

%% Problem 1.d 
% Plot
figure
hold on
scatter(xData,yAvg)
plot(xData, fit(pmin,xData), 'k')
plot(xData, fit([AVal_ols, pVal_ols],xData),'--k')
plot(xData, fit([AVal_wls, pVal_wls],xData),':k')
errorbar(xData,yAvg, ySEM, 'ok')
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([0,5]);
xlabel('X')
ylabel('Y')
title('y = Ax^p')
legend({'yAvg',sprintf('\\chi^2'),'OLS', 'WLS'})

%% End of script
fclose(fileName);
     