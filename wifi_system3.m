clear all;
% close all;

if exist('nameAndIdx.mat')
    delete('nameAndIdx.mat');
end
if exist('wifiStationAndIdx.mat')
    delete('wifiStationAndIdx.mat');
end

% path = './walter-basement-Feb26-part1/';
path = './walter-basement-Feb26-part2-data/';

path_wifi = [path, 'wifi/'];
files = dir(path_wifi);
%% The first two are '.' and '..'
file_names = {files(3:end).name}';

data = [];
id_names = {};
for i = 1:length(file_names)
    fid = fopen([path_wifi, file_names{i}]);
    token = textscan(fid, '%d64%d64%d64%f64%d64%s\n');
    
    data = [data; double((i-1)*ones(length(token{1}), 1)), wifiStation2idx(double(token{1})), double(token{2}), (10.^(double(token{2})./10)), double(token{3})*1e-6, double(token{4}), double(token{5}), double(name2idx(token{6}))];
    fclose(fid);
end

%% Remove Duplicate
[~, IA, ~] = unique(data(:,[2, 3, 5, 7, 8]), 'rows', 'stable');
%% create data matrix without duplicate
unique_data = data(IA, :);
for i=min(unique_data(:,1)):max(unique_data(:,1))
    if sum(unique_data(:,1)==i)==0
        continue;
    end
    unique_data(unique_data(:,1)==i & unique_data(:,5) < (max(unique_data(unique_data(:,1)==i,5)) - 1), :) = [];
end

for i = 1:max(unique_data(:,8))
    if length(unique(unique_data(unique_data(:,8)==i, 2))) <= 1
        unique_data(unique_data(:,8)==i, :) = [];
    end
end

f_min = min(unique_data(:,7));
unique_data(:,7) = unique_data(:,7) / f_min;

temp_data = [];
j = 1;
for i = 1:unique_data(end,1)
    curr_data = unique_data(unique_data(:,1)==i, [1:3,5:7]);

    if isempty(curr_data)
        continue;
    end
    
    curr_data(:,1) = j;
    aps = unique(curr_data(:,2));
    for k = 1:length(aps)
        ap_data = curr_data(curr_data(:,2)==aps(k),:);
        freqs = unique(ap_data(:,6));
        
        for l = 1:length(freqs)
            freq_data = ap_data(ap_data(:,6)==freqs(l), :);
            median_freq_data = median(freq_data,1);
            temp_data = [temp_data; median_freq_data];
        end
    end
    
    j = j+1;
end
unique_data = temp_data;

min_count = 10;
station_freq = find(hist(unique_data(:,2), max(unique_data(:,2))) <= min_count);
for i = 1:length(station_freq)
    unique_data(unique_data(:,2)==station_freq(i), :) = [];
end
[~,~,ic1] = unique(unique_data(:,1),'rows','stable');
C1 = (1:length(ic1))';
unique_data(:,1) = C1(ic1);
[C,~,ic] = unique(unique_data(:,2),'rows', 'stable');
unique_data(:,2) = ic;      %% unique_data(:,2) == C(ic)

wifi_timestamps = mean(unique_data(:,5)-unique_data(:,4))+unique_data(:,4);     %todo compare max

%% Time Sync
imu = dlmread([path, 'imu_interpolated.txt']);
imu = imu(:,1:2);
imu(:,2)= imu(:,2)*1e-3;
diff = imu(:,1)-imu(:,2);
mean_diff = mean(diff);
std_diff = std(diff);
mean_diff = mean(diff(diff < (mean_diff+2*std_diff) & diff > (mean_diff - 2*std_diff)));

wifi_timestamps = wifi_timestamps + mean_diff;

%% Process data
g_ba_ci = BA_interp_pos( path );

fid = fopen([path, 'timestamps.txt']);
token = textscan(fid, '%s%f64\n');
img_timestamps = token{2};

N = unique_data(end,1);
Ns = max(unique_data(:,2));
init_ba_indx = 5;

A = [];
b = [];
xyzfr = {};
mean_zs = [];

Pj = -30;
gamma = 2;

for i = 1:Ns
    curr_data = unique_data(unique_data(:,2)==i, [3, 1, 6]);
    d_sq = (10.^((Pj - curr_data(:,1))/5/gamma)) ./ (curr_data(:,3).^2);

    if isempty(curr_data)
        continue;
    end
    
    timestamps = wifi_timestamps(unique_data(:,2)==i);
    xi = zeros(length(timestamps), 1);
    yi = zeros(length(timestamps), 1);
    z_loc = zeros(length(timestamps), 1);
    
    for j = 1:length(timestamps)
        [mx,inds] = sort(abs(img_timestamps - timestamps(j)));
        inds = inds(1:2);
        times = img_timestamps(inds);
        if max(inds)-init_ba_indx+1 > length(g_ba_ci)
            continue;
        end
        if (times(2)-timestamps(j))*(times(2)-timestamps(j)) > 0
            c1_p_ci = g_ba_ci(:,inds(1)-init_ba_indx+1);
        else
            c1_p_ci = (g_ba_ci(:,inds(1)-init_ba_indx+1) * abs(times(2)-timestamps(j)) + g_ba_ci(:,inds(2)-init_ba_indx+1) * abs(times(1)-timestamps(j))) / abs(times(2)-times(1));
        end
        xi(j) = c1_p_ci(1);
        yi(j) = c1_p_ci(2);
        z_loc(j) = c1_p_ci(3);
    end

    xyzfr{i} = {[xi, yi, z_loc, curr_data(:,3), curr_data(:,2), curr_data(:,1)]};
    
    x_loc_sq = xi .* xi;
    y_loc_sq = yi .* yi;
    z_loc_sq = z_loc .* z_loc;

    a_vec = (x_loc_sq - mean(x_loc_sq)) + (y_loc_sq - mean(y_loc_sq)) + (z_loc_sq - mean(z_loc_sq));
    b_vec = 2 * (xi - mean(xi));
    c_vec = 2 * (yi - mean(yi));
    d_vec = 2 * (z_loc - mean(z_loc));
    e_vec = d_sq - mean(d_sq);
    mean_zs = [mean_zs, mean(z_loc)];
    
    A_temp = [b_vec, c_vec, d_vec, e_vec];
    A = [A, zeros(size(A,1), size(A_temp, 2)); zeros(size(A_temp,1), size(A, 2)), A_temp];
    b = [b; a_vec];
end

x = (A'*A)\(A'*b);
ap_loc = reshape(x, 4, length(x)/4);
median_scale = sqrt(median(ap_loc(4,:)));
mean_z = mean(mean_zs);
if median_scale < 0
    median_scale = 1;
end
x(3:4:end) = x(3:4:end) * median_scale - mean_z;
x(4:4:end) = Pj;
x(length(x)+1) = gamma;

figure, axis equal, hold on
plot3(g_ba_ci(1,:), g_ba_ci(2,:), g_ba_ci(3,:),'r')
for j = 1:length(ap_loc)
    plot3(ap_loc(1,j), ap_loc(2,j), ap_loc(3,j)*median_scale, 'bx')
end

err = 1e6;
prev_err = 1e9;
prev_x = x;
while (err < prev_err) && (prev_err-err > 1e-3)
%     J = zeros(N, 4*Ns+1);
    J = [];
%     r = zeros(N,1);
    r = [];
    count = zeros(N,1);
    
    for j = 1:Ns
        xyzfr_curr = cell2mat(xyzfr{j});
        xi = xyzfr_curr(:,1);
        yi = xyzfr_curr(:,2);
        fjk = xyzfr_curr(:,4);
%         index = xyzfr_curr(:,5);
        rssi = xyzfr_curr(:,6);
        
        xj = x(4*j-3);
        yj = x(4*j-2);
        zj = x(4*j-1);
        Pj = x(4*j);
        gamma = x(end);
        dij = sqrt((xj-xi).*(xj-xi) + (yj-yi).*(yj-yi) + zj.*zj);
        
        Jx = 2 * (xj - xi) ./ dij;
        Jy = 2 * (yj - yi) ./ dij;
        Jz = 2 * zj ./ dij;
        dd = 10 .^ ((Pj - rssi) / 10 / gamma);
        JP = -log(10)/10.0/gamma ./ fjk .* dd;
        Jn = log(10)/10.0/gamma/gamma ./ fjk .* dd;
        
%         J(index, 4*j-3) = J(index, 4*j-3) + Jx;
%         J(index, 4*j-2) = J(index, 4*j-2) + Jy;
%         J(index, 4*j-1) = J(index, 4*j-1) + Jz;
%         J(index, 4*j) = J(index, 4*j) + JP;
%         J(index, end) = J(index, end) + Jn;
%         
%         r(index) = r(index) + dij - dd ./ fjk;
%         count(index) = count(index) + 1;
        J_temp = zeros(length(Jx), 4*Ns+1);
        J_temp(:,4*j-3) = Jx;
        J_temp(:,4*j-2) = Jy;
        J_temp(:,4*j-1) = Jz;
        J_temp(:,4*j) = JP;
        J_temp(:,end) = Jn;
        
        r_temp = dij - dd ./ fjk;
        J = [J; J_temp];
        r = [r; r_temp];
    end
    
%     r = r ./ count;
    dx = (J'*J)\(J'*r);
    
    prev_x = x;
    x = x + dx;
    prev_err = err;
    err = sum(abs(r)) / length(r)
    
%     ap_loc = reshape(x(1:end-1), 4, length(x(1:end-1))/4);
%     figure, axis equal, hold on
%     plot3(g_ba_ci(1,:), g_ba_ci(2,:), zeros(1,length(g_ba_ci)),'r')
%     for j = 1:length(ap_loc)
%         plot3(ap_loc(1,j), ap_loc(2,j), ap_loc(3,j), 'bx')
%     end
%     keyboard
end

if err > prev_err
    x = prev_x;
end

ap_loc = reshape(x(1:end-1), 4, length(x(1:end-1))/4);
figure, axis equal, hold on
plot3(g_ba_ci(1,:), g_ba_ci(2,:), zeros(1,length(g_ba_ci)),'r')
for j = 1:length(ap_loc)
    plot3(ap_loc(1,j), ap_loc(2,j), ap_loc(3,j), 'bx')
end
