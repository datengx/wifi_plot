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
[~, IA, ~] = unique(data(:,[2, 5]), 'rows', 'stable');
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

min_count = 4;
station_freq = find(hist(unique_data(:,2), max(unique_data(:,2))) <= min_count);
for i = 1:length(station_freq)
    unique_data(unique_data(:,2)==station_freq(i), :) = [];
end
[~,~,ic1] = unique(unique_data(:,1),'rows','stable');
C1 = (1:length(ic1))';
unique_data(:,1) = C1(ic1);
[C,~,ic] = unique(unique_data(:,2),'rows', 'stable');
unique_data(:,2) = ic;      %% unique_data(:,2) == C(ic)

wifi_timestamps = mean(unique_data(:,6)-unique_data(:,5))+unique_data(:,5);     %todo compare max
% plot(diff(unique_data(:,6)),'b-o'),hold on, plot(diff(unique_data(:,5)), 'r-x')

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

X = zeros(3, Ns);
A = [];
b = [];
xyd = {};

for i = 1:Ns
    curr_data = unique_data(unique_data(:,2)==i, [3, 7, 6, 1]);
    d_sq = (((75.0 / pi) ./ curr_data(:,2)).^2) .* (10.^-(curr_data(:,1)/10));
%     d_sq = 10.^-((curr_data(:,1)+20*log(curr_data(:,2)))/10);

%     curr_data = unique_data(unique_data(:,2)==i, [4,7, 6]);
%     d_sq = 1 ./ (curr_data(:,1));

    if isempty(curr_data)
        continue;
    end
    
    timestamps = wifi_timestamps(unique_data(:,2)==i);
    x_loc = zeros(length(timestamps), 1);
    y_loc = zeros(length(timestamps), 1);
    
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
    end

    xyd{i} = {[x_loc, y_loc, d_sq, curr_data(:,4)]};
    
    x_loc_sq = x_loc .* x_loc;
    y_loc_sq = y_loc .* y_loc;

    a_vec = (x_loc_sq - mean(x_loc_sq)) + (y_loc_sq - mean(y_loc_sq));
    b_vec = 2 * (x_loc - mean(x_loc));
    c_vec = 2 * (y_loc - mean(y_loc));
    d_vec = d_sq - mean(d_sq);
    
    A_temp = [b_vec, c_vec, d_vec];
    A = [A, zeros(size(A,1), size(A_temp, 2)); zeros(size(A_temp,1), size(A, 2)), A_temp];
    b = [b; a_vec];
end

x = (A'*A)\(A'*b);
for i = 3:3:length(x)
    x(i) = real(10 * log10(x(i)));
end
ap_loc = reshape(x, 3, length(x)/3);

figure, axis equal, hold on
plot(g_ba_ci(1,:), g_ba_ci(2,:),'r')
for j = 1:length(ap_loc)
    plot(ap_loc(1,j), ap_loc(2,j), 'bx')
end

err = 1e6;
while err > 1e-3
    J = zeros(N, 3*Ns);     %todo
    r = zeros(N,1);
    count = zeros(N,1);
    
    for i = 1:Ns
        xyd_curr = cell2mat(xyd{i});
        x_loc = xyd_curr(:,1);
        y_loc = xyd_curr(:,2);
        d_sq = xyd_curr(:,3);
        index = xyd_curr(:,4);
        
        xj = x(3*i-2);
        yj = x(3*i-1);
        wj = x(3*i);
        
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