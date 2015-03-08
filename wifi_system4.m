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

% f_min = min(unique_data(:,7));
% unique_data(:,7) = unique_data(:,7) / f_min;

unique_data(unique_data(:,3)<=-85,:) = [];
temp_data = [];
n = [];
% j = 1;
for i = 1:unique_data(end,1)
    curr_data = unique_data(unique_data(:,1)==i, [1:3,5:7]);
    
    if isempty(curr_data)
        continue;
    end
    
    %     curr_data(:,1) = j;
    aps = unique(curr_data(:,2));
    for k = 1:length(aps)
        ap_data = curr_data(curr_data(:,2)==aps(k),:);
        freqs = unique(ap_data(:,6));
        
        for l = 1:length(freqs)
            freq_data = ap_data(ap_data(:,6)==freqs(l), :);
            median_freq_data = median(freq_data,1);
            temp_data = [temp_data; median_freq_data];
        end
%         if length(freqs)==2
%             r_curr = temp_data(end-1:end,3);
%             f_curr = temp_data(end-1:end,6);
%             n_curr = (r_curr(1)-r_curr(2)) / log10(f_curr(2)/f_curr(1));
%             n = [n; n_curr/10];
%         end
n = [n; length(freqs)];
    end
    
    %     j = j+1;
end
unique_data = temp_data;

min_count = 10;
station_freq = find(hist(unique_data(:,2), max(unique_data(:,2))) <= min_count);
for i = 1:length(station_freq)
    unique_data(unique_data(:,2)==station_freq(i), :) = [];
end
% [~,~,ic1] = unique(unique_data(:,1),'rows','stable');
% C1 = (1:length(ic1))';
% unique_data(:,1) = C1(ic1);
% [C,~,ic] = unique(unique_data(:,2),'rows', 'stable');
% unique_data(:,2) = ic;      %% unique_data(:,2) == C(ic)

wifi_timestamps = mean(unique_data(:,5)-unique_data(:,4))+unique_data(:,4);     %todo compare max

%% Time Sync
imu = dlmread([path, 'imu_interpolated.txt']);
imu = imu(:,1:2);
imu(:,2)= imu(:,2)*1e-3;
diff_imu = imu(:,1)-imu(:,2);
mean_diff = mean(diff_imu);
std_diff = std(diff_imu);
mean_diff = mean(diff_imu(diff_imu < (mean_diff+2*std_diff) & diff_imu > (mean_diff - 2*std_diff)));

wifi_timestamps = wifi_timestamps + mean_diff;

%% Process data
g_ba_ci = BA_interp_pos( path );

fid = fopen([path, 'timestamps.txt']);
token = textscan(fid, '%s%f64\n');
img_timestamps = token{2};

init_ba_indx = 5;

for j = 1:size(wifi_timestamps,1)
    [mx,inds] = sort(abs(img_timestamps - wifi_timestamps(j)));
    inds = inds(1:2);
    times = img_timestamps(inds);
    if max(inds)-init_ba_indx+1 > length(g_ba_ci)
        continue;
    end
    if (times(2)-wifi_timestamps(j))*(times(2)-wifi_timestamps(j)) > 0
        c1_p_ci = g_ba_ci(:,inds(1)-init_ba_indx+1);
    else
        c1_p_ci = (g_ba_ci(:,inds(1)-init_ba_indx+1) * abs(times(2)-wifi_timestamps(j)) + g_ba_ci(:,inds(2)-init_ba_indx+1) * abs(times(1)-wifi_timestamps(j))) / abs(times(2)-times(1));
    end
    
    xyz(j,:) = c1_p_ci';
end

unique_data(:,4:5) = [];
unique_data = [unique_data, wifi_timestamps, xyz];

basement_ap_ids = [5,6,10,14,20,21,22];
temp_data = [];
for i = 1:length(basement_ap_ids)
    temp_data = [temp_data; unique_data(unique_data(:,2)==basement_ap_ids(i), :)];
end
unique_data = temp_data;

%% Position of WIFI station in the basement
idAndPosition = [5, -2.66,    11.8994, 1.2192;...
                 6,  13.106, -19.62,   1.2192;...
                10, -15.27,  -13.22,   1.1684;...
                14, -15.6179,-32.029,  1.27;...
                20, -24.054, -2.66,    1.1684;...
                21, -47.34,  -11.22,   1.4732;...
                22, -45.58,  -27.72,   1.4982];
            
%% Least square to solve n and A
[N,M] = size(unique_data);
unique_data = [unique_data,zeros(N,3)];
% z_s_avg = unique_data(:,8)/N;
for k = 1:N
    unique_data(k,9:11) = idAndPosition(...
                         idAndPosition(:,1) == unique_data(k,2),2:4);
end

R = unique_data(:,3);
D = 10*log10( sqrt((unique_data(:,6)-unique_data(:,9)).^2 ...
                 + (unique_data(:,7)-unique_data(:,10)).^2 ...
                 + (unique_data(:,8)-unique_data(:,11)).^2)...
               .* unique_data(:,4));

n_A = [D,-ones(N,1)] \ (-R)

% figure, axis equal, hold on
% plot3(g_ba_ci(1,:), g_ba_ci(2,:), g_ba_ci(3,:),'r')
% ap_ids = unique(unique_data(:,2));
% colors = jet(length(ap_ids));
% for i = 7
%     j = ap_ids(i);
%     curr_data = unique_data(unique_data(:,2)==j,:);
%     freqs = unique(curr_data(:,4));
%     for k = 1:size(curr_data,1)
%         if curr_data(k,4) == freqs(1)
%             continue;
%         end
%         d = 10.^(-curr_data(k,3)/19) *.001;
%         circle(curr_data(k,6), curr_data(k,7), d);
%     end
%     scatter3(curr_data(:,6), curr_data(:,7), curr_data(:,8)+2*i, 'MarkerEdgeColor', colors(i,:), 'Marker', 'x')
% 
% %     figure, hold on, axis equal
% %     for l = 1:length(freqs)
% %         scatter3(curr_data(curr_data(:,4)==freqs(l),6),curr_data(curr_data(:,4)==freqs(l),7),curr_data(curr_data(:,4)==freqs(l),3),'MarkerEdgeColor', colors(l*floor(size(colors,1)/2),:), 'Marker', '.')
% %     end
% end

