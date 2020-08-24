clear all;
close all;
clc;
format long g;

%% Inputs
prn = 11;

preciseFile = 'igs20872.sp3';
brdc24File = 'brdc24hr.txt';
brdcbestFile = 'brdcbest.txt';

almanacFile = 'Almanac-SatPos-PRN11.txt';
urPredFile = 'igu20872_00.sp3';
urObsFile = 'igu20873_00.sp3';

% Roof N1 pillar coordinates of receiver
N1 = [-1641900.798; -3664874.335; 4939969.355];
lat = 51.07942678;
lon = -114.13286260;

%% Read Files
refsat=readSP3(preciseFile);
epochs = unique(refsat(:, 3));
PRN = unique(refsat(:, 1));
refsat(:, 4:6) = refsat(:, 4:6)*1000;
%% ------------------ PART 1------------------
% Sky plot of satellite
   %calculate the vector from receiver to satellites in xyz
   %transform the vector to ENU
   %calculate the azimuth and elevation
   %plot the sky plot
   %calculate the DOP
   
%% ------------------ TASK 1.1 ------------------

%receiver to sat vector in enu 
vector_rs = zeros(length(refsat),3);
vector_rs_ENU = zeros(length(refsat),3);
for i = 1:length(refsat)
    vector_rs(i,:) = [refsat(i,4) - N1(1), refsat(i,5) - N1(2), refsat(i,6) - N1(3)]; 
    currentVec = vector_rs(i,:)';
    vector_rs_ENU(i,:) = xyz2enu(lat, lon, currentVec);
end
   
% Azimuth and elevation angles
elevation = zeros(length(refsat), 1);
az = zeros(length(refsat), 1);
for i = 1:length(refsat)
    elevation(i) = asind(vector_rs_ENU(i,3)/(sqrt(vector_rs_ENU(i,1)^2 + vector_rs_ENU(i,2)^2 + vector_rs_ENU(i,3)^2))); % elevation angle [dec deg]
    az(i)   = atan2d(vector_rs_ENU(i,1),vector_rs_ENU(i,2)); % azimuth [dec deg]
    if elevation(i) <= 0
        elevation(i) = NaN;
    end
end

% Store angles in the RefSat matrix
refsat(:, 8) = elevation';
refsat(:, 9) = az';

count = 1; 
elev = zeros(length(epochs), length(PRN));
azimuth = zeros(length(epochs), length(PRN));

for i = 1:length(epochs)
    for j = 1:length(PRN)
        elev(i, j) = refsat(count, 8);
        azimuth(i, j) = refsat(count, 9);
        count = count + 1; 
    end
end

% Sky plot
figure;
for i = 1:length(PRN)
    skyplot(azimuth(:,i), elev(:,i), '-');
    hold on; 
end
title ({'Visible satellites from N1 pillar in a  24-hr period, January 07, 2020'}); 
hold off;
   

%% ------------------ TASK 1.2 ------------------
count = 1; 

for j = 1:length(epochs)
    stats(j).satNum(length(PRN), 1) = zeros();  % numSat = no. of visible satellites (32x1)
    stats(j).xyzK(length(PRN),3)   = zeros();   % xyzK = satellite coords (32x3)
    for i = 1:length(PRN)
        if ~isnan(elev(j,i))
            stats(j).satNum(i) = i;
            stats(j).xyzK(i,:) = refsat(count, 4:6);
        end
        count = count +1; 
    end
    % Calculate the design matrix A of N1 at each epoch 
    stats(j).A = design_matrix(stats(j).xyzK(:,:), N1(1:3,1));
    % Cofactor Matrix Qxyz   
    stats(j).Qxyz = ((stats(j).A)'*stats(j).A)^-1;
    % Number of visible satellites
    vis_sats(j) = size(stats(j).A,1); 
    % Q convertion to enu
    stats(j).Qenu = Qxyz2Qenu(lat, lon, stats(j).Qxyz);
    
    % DOPs
    EDOP(j,1) = sqrt(stats(j).Qenu(1,1));
    NDOP(j,1) = sqrt(stats(j).Qenu(2,2));
    HDOP(j,1)  = sqrt(stats(j).Qenu(1,1) + stats(j).Qenu(2,2));
    VDOP(j,1)  = sqrt(stats(j).Qenu(3,3));

end
   
% Plots
t = 1:length(epochs);
t = t';

% Plotting the DOPs
figure
subplot(3,1,[1 2])
plot(t, VDOP(:,1));
hold on; 
plot(t, HDOP(:,1));
plot(t, NDOP(:,1));
plot(t, EDOP(:,1));
legend ('VDOP','HDOP','NDOP','EDOP');
title({'DOP values and number of visible satellites versus epoch'});
xlabel('Epoch (900 second intervals)');
ylabel('DOP');
ylim ([0,2]);
grid on;
hold off;

subplot(313)
plot (t, vis_sats); 
xlabel('Epoch (900 second intervals)');
ylabel('No. of Satellites');
ylim ([8,16]);
grid on;

%% ------------------ TASK 3 ------------------
%% ------------------ TASK 3.1 ------------------
file_best = load(brdcbestFile);
file_24 = load(brdc24File);

%% ------------------ TASK 3.2 ------------------
almanac = load(almanacFile);

%% ------------------ TASK 3.3 ------------------
file_predicted = readSP3(urPredFile);
file_observed = readSP3(urObsFile);
idx_predicted = find(file_predicted(:,1) == prn);
idx_observed = find(file_observed(:,1) == prn);
xyz_predicted = file_predicted(idx_predicted,:);
xyz_observed = file_observed(idx_observed,:);

% Convert to km
xyz_predicted(:, 4:6) = xyz_predicted(:, 4:6)*1000;
xyz_observed(:, 4:6) = xyz_observed(:, 4:6)*1000;

%% ------------------ TASK 3.4 ------------------
file_precise = readSP3(preciseFile);
idx_precise = find(file_precise(:,1) == prn);
xyz_precise = file_precise(idx_precise,:);
xyz_precise(:, 4:6) = xyz_precise(:, 4:6)*1000;

%% ------------------ TASK 3.5 ------------------

error_24 = zeros(length(file_24),3);
error_best = zeros(length(file_best),3);
error_almanac = zeros(length(almanac),3);
error_predicted = zeros(length(file_best),3);
error_observed = zeros(length(file_best),3);
for i = 1:length(xyz_precise)
    for j = 1:length(file_24)
        if xyz_precise(i,3) == file_24(j,1)
            error_24(i, 1) = file_24(j,2)-xyz_precise(i,4);
            error_24(i, 2) = file_24(j,3)-xyz_precise(i,5);
            error_24(i, 3) = file_24(j,4)-xyz_precise(i,6);
         end
    end
    for j = 1:length(file_best)
        if xyz_precise(i,3) == file_best(j,1)
            error_best(i, 1) = file_best(j,2)-xyz_precise(i,4);
            error_best(i, 2) = file_best(j,3)-xyz_precise(i,5);
            error_best(i, 3) = file_best(j,4)-xyz_precise(i,6);
        end
    end
    for j = 1:length(almanac)
        if xyz_precise(i,3) == almanac(j,1)
            error_almanac(i, 1) = almanac(j,2)-xyz_precise(i,4);
            error_almanac(i, 2) = almanac(j,3)-xyz_precise(i,5);
            error_almanac(i, 3) = almanac(j,4)-xyz_precise(i,6);
        end
    end
    for j = 1:length(xyz_predicted)
        if xyz_precise(i,3) == xyz_predicted(j,3)
            error_predicted(i, 1) = xyz_predicted(j,4)-xyz_precise(i,4);
            error_predicted(i, 2) = xyz_predicted(j,5)-xyz_precise(i,5);
            error_predicted(i, 3) = xyz_predicted(j,6)-xyz_precise(i,6);
        end
    end
    for j = 1:length(xyz_observed)
        if xyz_precise(i,3) == xyz_observed(j,3)
            error_observed(i, 1) = xyz_observed(j,4)-xyz_precise(i,4);
            error_observed(i, 2) = xyz_observed(j,5)-xyz_precise(i,5);
            error_observed(i, 3) = xyz_observed(j,6)-xyz_precise(i,6);
        end
    end
end

%% ------------------ TASK 3.6 ------------------

figure 
plot(t, error_24(:,1));
hold on
plot(t, error_24(:,2));
plot(t, error_24(:,3));
hold off
grid on
title ({'PRN 11 position error with broadcast ephemeris -','Precision compared to 24-hour solution data'}); 
legend ('Error in x', 'Error in y', 'Error in z'); 
xlabel('Epoch (900-second intervals)');
ylabel('Error (metres)');

figure 
plot(t, error_best(:,1));
hold on
plot(t, error_best(:,2));
plot(t, error_best(:,3));
hold off
grid on
title ({'PRN 11 position error with broadcast ephemeris -','Precision compared to best solution data'}); 
legend ('Error in x', 'Error in y', 'Error in z'); 
xlabel('Epoch (900-second intervals)');
ylabel('Error (metres)');

figure 
plot(t, error_almanac(:,1));
hold on
plot(t, error_almanac(:,2));
plot(t, error_almanac(:,3));
hold off
grid on
title ({'PRN 11 position error with broadcast ephemeris -','Precision compared to Almanac data'}); 
legend ('Error in x', 'Error in y', 'Error in z'); 
xlabel('Epoch (900-second intervals)');
ylabel('Error (metres)');

figure 
plot(t, error_predicted(:,1));
hold on
plot(t, error_predicted(:,2));
plot(t, error_predicted(:,3));
hold off
grid on
title ({'PRN 11 position error with broadcast ephemeris -','Precision compared to IGS predicted data'}); 
legend ('Error in x', 'Error in y', 'Error in z'); 
xlabel('Epoch (900-second intervals)');
ylabel('Error (metres)');

figure 
plot(t, error_observed(:,1));
hold on
plot(t, error_observed(:,2));
plot(t, error_observed(:,3));
hold off
grid on
title ({'PRN 11 position error with broadcast ephemeris -','Precision compared to IGS observed data'}); 
legend ('Error in x', 'Error in y', 'Error in z');
xlabel('Epoch (900-second intervals)');
ylabel('Error (metres)');

%% ------------------ TASK 4 ------------------
%------------------Statistics------------------

% Maximum Absolute Error
eMax(1,1:3) = [max(abs(error_best(:,1))),max(abs(error_best(:,2))),max(abs(error_best(:,3)))]
eMax(2,1:3) = [max(abs(error_24(:,1))),max(abs(error_24(:,2))),max(abs(error_24(:,3)))]
eMax(3,1:3) = [max(abs(error_almanac(:,1))),max(abs(error_almanac(:,2))),max(abs(error_almanac(:,3)))]
eMax(4,1:3) = [max(abs(error_predicted(:,1))),max(abs(error_predicted(:,2))),max(abs(error_predicted(:,3)))]
eMax(5,1:3) = [max(abs(error_observed(:,1))),max(abs(error_observed(:,2))),max(abs(error_observed(:,3)))]

% Mean Error
eMean(1,1:3) = [mean(error_best(:,1)), mean(error_best(:,2)), mean(error_best(:,3))]
eMean(2,1:3) = [mean(error_24(:,1)), mean(error_24(:,2)), mean(error_24(:,3))]
eMean(3,1:3) = [mean(error_almanac(:,1)), mean(error_almanac(:,2)), mean(error_almanac(:,3))]
eMean(4,1:3) = [mean(error_predicted(:,1)), mean(error_predicted(:,2)), mean(error_predicted(:,3))]
eMean(5,1:3) = [mean(error_observed(:,1)), mean(error_observed(:,2)), mean(error_observed(:,3))]

% Standard Deviation 
eStd(1,1:3) = [std(error_best(:,1)), std(error_best(:,2)), std(error_best(:,3))]
eStd(2,1:3) = [std(error_24(:,1)), std(error_24(:,2)), std(error_24(:,3))]
eStd(3,1:3) = [std(error_almanac(:,1)), std(error_almanac(:,2)), std(error_almanac(:,3))]
eStd(4,1:3) = [std(error_predicted(:,1)), std(error_predicted(:,2)), std(error_predicted(:,3))]
eStd(5,1:3) = [std(error_observed(:,1)), std(error_observed(:,2)), std(error_observed(:,3))]

% RMSE
eRMS(1,1:3) = [rms(error_best(:,1)), rms(error_best(:,2)), rms(error_best(:,3))]
eRMS(2,1:3) = [rms(error_24(:,1)), rms(error_24(:,2)), rms(error_24(:,3))]
eRMS(3,1:3) = [rms(error_almanac(:,1)), rms(error_almanac(:,2)), rms(error_almanac(:,3))]
eRMS(4,1:3) = [rms(error_predicted(:,1)), rms(error_predicted(:,2)), rms(error_predicted(:,3))]
eRMS(5,1:3) = [rms(error_observed(:,1)), rms(error_observed(:,2)), rms(error_observed(:,3))]

%% ------------------ TASK 5 ------------------

avgHDOP = mean(HDOP)
avgVDOP = mean(VDOP)

%% Functions
function [ enu ] = xyz2enu (lat, lon, vector)
% Rotation matrix, R
R = [-sind(lon),cosd(lon),0;
 -sind(lat)*cosd(lon),-sind(lat)*sind(lon),cosd(lat);
 cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat)];

% xyz to enu
enu = R*vector;
enu = enu';
end

function [ Q ] = Qxyz2Qenu( lat, lon, Qx )
% Rotation matrix, R
R = [-sind(lon),cosd(lon),0,0;
 -sind(lat)*cosd(lon),-sind(lat)*sind(lon),cosd(lat), 0;
 cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat), 0;
 0,0,0,1];

% xyz to enu
Q = R*Qx*(R');
end

function [ A ] = design_matrix( xyz_sat, xyzR )
satnum = 1; 
for i = 1:length(xyz_sat)
    if xyz_sat(i, 1) == 0
        continue
    else
        % Numerator
        dx = xyzR(1, 1) - xyz_sat(i, 1);
        dy = xyzR(2, 1) - xyz_sat(i, 2);
        dz = xyzR(3, 1) - xyz_sat(i, 3);
        % geometric distance from reciever to satellite
        d = sqrt(dx^2+dy^2+dz^2);
        % Design matrix A
        A(satnum, :) = [dx/d, dy/d, dz/d, -1];
        satnum = satnum + 1; 
    end
end 

end
