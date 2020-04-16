%% Data analysis for Major Experiment
close all;
load data;
showPlot = true;
showContour = false;

%% Constants for measurement uncertainties
TEMP_READING_UNCERTAINTY = 0.5; % Celsius 
TEMP_VOLTAGE_UNCERTAINTY = 9.1e-3; % Volts
OPTICAL_READING_UNCERTAINTY = 5e-3; % mW
OPTICAL_VOLTAGE_UNCERTAINTY = 9.1e-3; % Volts
%% Organize the data for binning (only temp)
tempX = tempData(:,1);
tempY = tempData(:,2);

% Loop through the data, organize all y measurements per each x
string = "";
i = 2;
while i < length(tempY)
    string = strcat(string,num2str(tempY(i-1)), "\t");
    while i <= 203 && tempY(i-1) == tempY(i)
         string = strcat(string, num2str(tempX(i-1)),"\t");

         i = i + 1; 
    end
    string = strcat(string,"\n");
    i = i + 1;
end

% Compute average x for each y value
tempY_avg = zeros(length(string), 1);
tempX_avg = zeros(length(string), 1);
tempX_std = zeros(length(string), 1);
tempX_nVals = zeros(length(string), 1);
string = split(string,"\n");
for i=1:length(string)-1
   tempString = split(string(i),"\t");
   tempY_avg(i) = str2double(tempString(1));
   % Find the average x values
   xSum = 0;
   xLen = length(tempString) -2;
   for j=2:length(tempString)-1
       xSum = xSum + str2double(tempString(j)); 
   end
   xAvg = xSum / (xLen);
   % Find the standard deviation of x values
   dMean = 0;
   for j=2:length(tempString)-1
      dMean = dMean + power((str2double(tempString(j)) - xAvg),2); 
   end
   xStd = sqrt(dMean/xLen);
   % Add all statistics to their respective arrays
   tempX_avg(i) = xAvg;
   tempX_std(i) = xStd + TEMP_VOLTAGE_UNCERTAINTY;
   tempX_nVals(i) = xLen;
end

%% Chi^2 Analysis: Temperature TODO: Uncertainty in each parameter.

% Find SEM for x values (voltage)
tempYStd = zeros(length(tempY_avg),1) + TEMP_READING_UNCERTAINTY;
tempYSEM = tempYStd' ./ sqrt(tempX_nVals);
% tempFit for thermistor
% Equation: T(V) = Beta / [ ln(Rs/Rinf) * (V+/V -1) ] - 273.15
tempFit = @(p, tempX_avg) p(1) ./ (log(499 ./ p(2) .* (5 ./ ...
                tempX_avg - 1))) - 273.15;
chi2_temp = @(p) sum( ((tempY_avg - tempFit(p,tempX_avg)) ./ tempYSEM).^2);

pmin_temp = fminsearch(chi2_temp, [1000 0]);
chi2_temp_min = chi2_temp(pmin_temp);
chi2_temp_reduced = chi2_temp_min / (length(tempX_avg) - 2);
fprintf('Beta = %f, R_inf = %f \n', pmin_temp(1), pmin_temp(2));
fprintf('Chi^2 Reduced: %f\n', chi2_temp_reduced);

%% Chi^2 Uncertainty: Temperature
% Find uncertainty in Beta and R_inf using Chi^2 + 1
Beta = pmin_temp(1);
R_inf = pmin_temp(2);
rB = 1;
rT = .001;
BRange = linspace(Beta - rB, Beta + rB, 100);
TRange = linspace(R_inf - rT, R_inf + rT, 100);

% 1. Generate a grid to search parameter pairs for min Chi^2 + 1
[BSpace, TSpace] = ndgrid(BRange, TRange); 
chiSquareGrid = zeros(size(BSpace));

% 2. In the grid, compute Chi^2 for every parameter combination
for i = 1:size(BSpace,1)
    for j = 1:size(BSpace,2)
        chiSquareGrid(i,j) = chi2_temp([BSpace(i,j) TSpace(i,j)]);
    end
end    

% 3. Find the indices where chi^2 <= chi2_min + 1
[indexB, indexT] = find(chiSquareGrid <= chi2_temp_min + 1);

% 4. Find the min and max of these parameters
dBMin = min(BRange(indexB));
dBMax = max(BRange(indexB));
dTMin = min(TRange(indexT));
dTMax = max(TRange(indexT));

% 5. Compute the uncertainties!
dB = (dBMax - dBMin)/2;
dT = (dTMax - dTMin)/2;
fprintf('Beta Uncertainty: %f \t R_inf Uncertainty: %f \n', dB, dT);
if showContour
    figure
    hold on
    contourf(BSpace, TSpace, real(chiSquareGrid), 50);
    scatter(BRange(indexB), TRange(indexT), 'mo', 'MarkerFaceColor', 'Magenta');
    plot(pmin_temp(1), pmin_temp(2), 'kx', 'MarkerSize', 15);
    xlabel('\\Beta'); ylabel('R_inf'); ch = colorbar; ylabel(ch, ... 
        sprintf('\\chi^2'));
    title(sprintf('\\chi^2 Parameter Space'));
end
%% Organize Photodiode data for chi^2 fitting. 
photoX = opticalData(:,1);
photoY = opticalData(:,2);

photoYStd = zeros(length(photoY), 1) + OPTICAL_READING_UNCERTAINTY;
photoYSEM = photoYStd;

%% Chi^2 Analysis: Photodiode 

% Fit for Photodiode
% Equation: Ax^p
photoFit = @(p, photoX) p(1) .* (photoX) .^p(2);
chi2_photo = @(p) sum( ((photoY - photoFit(p,photoX)) ./ ...
                        1).^2);
pmin_photo = fminsearch(chi2_photo, [0,18]);
chi2_photo_min = chi2_photo(pmin_photo);
chi2_photo_reduced = chi2_photo_min / (length(photoX) - 2);
fprintf('A = %e, P = %f \n', pmin_photo(1), pmin_photo(2));
fprintf('Chi^2 Reduced: %f\n', chi2_photo_reduced);

%% Chi^2 Uncertainty: Photodiode
% Find uncertainty in A and p using Chi^2 + 1
A = pmin_photo(1);
p = pmin_photo(2);
rA = 1e-11;
rP = 1;
ARange = linspace(A - rA, A + rA, 100);
pRange = linspace(p - rP, p + rP, 100);

% 1. Generate a grid to search parameter pairs for min Chi^2 + 1
[ASpace, pSpace] = ndgrid(ARange, pRange); 
chiSquareGrid = zeros(size(ASpace));

% 2. In the grid, compute Chi^2 for every parameter combination
for i = 1:size(ASpace,1)
    for j = 1:size(ASpace,2)
        chiSquareGrid(i,j) = chi2_photo([ASpace(i,j) pSpace(i,j)]);
    end
end    

% 3. Find the indices where chi^2 <= chi2_min + 1
[indexA, indexP] = find(chiSquareGrid <= chi2_photo_min + 1);

% 4. Find the min and max of these parameters
dAMin = min(ARange(indexA));
dAMax = max(ARange(indexA));
dpMin = min(pRange(indexP));
dpMax = max(pRange(indexP));

% 5. Compute the uncertainties!
dA = (dAMax - dAMin)/2;
dP = (dpMax - dpMin)/2;
fprintf('A Uncertainty: %e \t P Uncertainty: %f \n', dA, dP);
if showContour
    figure
    hold on
    contourf(ASpace, pSpace, chiSquareGrid, 50);
    scatter(ARange(indexA), pRange(indexP), 'mo', 'MarkerFaceColor', 'Magenta');
    plot(pmin_photo(1), pmin_photo(2), 'kx', 'MarkerSize', 15);
    xlabel('A'); ylabel('P'); ch = colorbar; ylabel(ch, sprintf('\\chi^2'));
    title(sprintf('\\chi^2 Parameter Space'));
end

%% Make the x error bars
xErrorBar_temp = zeros(length(tempX_avg),1) +  TEMP_VOLTAGE_UNCERTAINTY;
xErrorBar_photo = zeros(length(photoX),1) + OPTICAL_VOLTAGE_UNCERTAINTY;
                    
%% Plot
if showPlot    
    figure
    hold on
    xlabel('Voltage (V)')
    ylabel('Temperature (C)')
    title('Thermistor Calibration')
    plot(tempX_avg,tempFit(pmin_temp,tempX_avg),'k')
    errorbar(tempX_avg,tempY_avg, tempYStd, tempYStd, xErrorBar_temp,...
        xErrorBar_temp,'ok')
    legend({sprintf(['T(V) = \\beta / [ ln(499/R_{inf}(5/V -1)) ] -' , ... 
        '273.15\n \\beta = %f %s %f ^oC\n R_{inf} = %f %s %f \\Omega'],...
        pmin_temp(1), char(177), dB, pmin_temp(2), char(177), dT)})

    hold off
    figure
    hold on
    title('Photoresistor Calibration')
    plot(photoX,photoFit(pmin_photo,photoX), 'k')
    errorbar(photoX, photoY, photoYSEM, photoYSEM, xErrorBar_photo,...
        xErrorBar_photo, 'ok')
    xlabel('Voltage (V)')
    ylabel('Irradiance (mW/cm^2)')
    legend({sprintf('y = Ax^p \n A = %e %s %e \n p = %f %s %f', ...
        pmin_photo(1), char(177), dA, pmin_photo(2), char(177), dP)})
    hold off
    
    figure 
    hold on
    title('Power vs. Voltage Output Comparison')
    plot(photoX, photoFit(pmin_photo, photoX) * .97^2, 'k')
    A = 884;
    B = -1.0700;
    R2 = 10000;
    vin = 5;
    vout = linspace(0.01,5);
    P = (R2/A)^(1/B) * (vin./vout - 1).^(1/B);
    plot(vout, P, '--k')
    xlabel('Voltage (V)')
    ylabel('Power (mW)')
    legend({'Our Fit', 'Fit from Tavares et al.'})
end



