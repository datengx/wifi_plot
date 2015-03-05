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
xyzf = {};
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
    x_loc = zeros(length(timestamps), 1);
    y_loc = zeros(length(timestamps), 1);
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
        x_loc(j) = c1_p_ci(1);
        y_loc(j) = c1_p_ci(2);
        z_loc(j) = c1_p_ci(3);
    end

    xyzf{i} = {[x_loc, y_loc, z_loc, curr_data(:,3), curr_data(:,2)]};
    
    x_loc_sq = x_loc .* x_loc;
    y_loc_sq = y_loc .* y_loc;
    z_loc_sq = z_loc .* z_loc;

    a_vec = (x_loc_sq - mean(x_loc_sq)) + (y_loc_sq - mean(y_loc_sq)) + (z_loc_sq - mean(z_loc_sq));
    b_vec = 2 * (x_loc - mean(x_loc));
    c_vec = 2 * (y_loc - mean(y_loc));
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

figure, axis equal, hold on
plot3(g_ba_ci(1,:), g_ba_ci(2,:), g_ba_ci(3,:),'r')
for j = 1:length(ap_loc)
    plot3(ap_loc(1,j), ap_loc(2,j), ap_loc(3,j)*median_scale, 'bx')
end

err = 1e6;
while err > 1e-3
    J = zeros(N, 4*Ns+1);
    r = zeros(N,1);
    count = zeros(N,1);
    
    for i = 1:Ns
        xyzf_curr = cell2mat(xyzf{i});
        x_loc = xyzf_curr(:,1);
        y_loc = xyzf_curr(:,2);
        f_ap = xyzf_curr(:,4);
        index = xyzf_curr(:,5);
        
        xj = x(4*i-3);
        yj = x(4*i-2);
        zj = x(4*i-1);
        
        Jx = -2 * (x_loc - xj);
        Jy = -2 * (y_loc - yj);
        Jw = -log(10)/10.0 * d_sq * (10^(wj/10));
        
        J(index, 3*i-2) = Jx;
        J(index, 3*i-1) = Jy;
        J(index, 3*i) = Jw;
        
        r(index) = r(index) + (x_loc-xj).^2 + (y_loc-yj).^2 - (10^(wj/10)) .* d_sq;
        count(index) = count(index) + 1;
    end
    
    r = r ./ count;
    dx = (J'*J)\(J'*r);
    
    x = x + dx;
    err = sum(abs(r))
    
    ap_loc = reshape(x, 3, length(x)/3);
    figure, axis equal, hold on
    plot(g_ba_ci(1,:), g_ba_ci(2,:),'r')
    for j = 1:length(ap_loc)
        plot(ap_loc(1,j), ap_loc(2,j), 'bx')
    end
end