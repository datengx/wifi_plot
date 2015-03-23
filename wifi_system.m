clear all;
% close all;

if exist('nameAndIdx.mat')
    delete('nameAndIdx.mat');
end
if exist('wifiStationAndIdx.mat')
    delete('wifiStationAndIdx.mat');
end

path = '/home/mars/desktop/Wifi/walter-basement-Feb26-part1/';

path_wifi = [path, 'wifi/'];
files = dir(path_wifi);
%% The first two are '.' and '..'
file_names = {files(3:end).name}';

data = [];
id_names = {};
for i = 1:length(file_names)
    fid = fopen([path_wifi, file_names{i}]);
    token = textscan(fid, '%d64%d64%d64%f64%d64%s\n');
    
    data = [data; double((i-1)*ones(length(token{1}), 1)), wifiStation2idx(double(token{1})), double(token{2}), (10.^(double(token{2})./10))*1e6, double(token{3})*1e-6, double(token{4}), double(token{5}), double(name2idx(token{6}))];
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

min_count = 4;
station_freq = find(hist(unique_data(:,2), max(unique_data(:,2))) <= min_count);
for i = 1:length(station_freq)
    unique_data(unique_data(:,2)==station_freq(i), :) = [];
end
% [~,~,ic1] = unique(unique_data(:,1),'rows','stable');
% C1 = (1:length(ic1))';
% unique_data(:,1) = C1(ic1);
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

N = unique_data(end,1) + 1;
Ns = max(unique_data(:,2));
init_ba_indx = 5;
test_data_num = 10;
dim = 2;

if dim == 2
    A = zeros(2*N, 2*Ns);
    b = zeros(2*N, 1);
    timestamps = zeros(N, 1);
    j = 1;
    
    for i = 1:N
      curr_data = unique_data(unique_data(:,1)==i-1, [2,3,7]); 
      curr_data(:,2)=10.^-((curr_data(:,2)+20*log(curr_data(:,3)))/20);
%         curr_data = unique_data(unique_data(:,1)==i-1, [2,4]);
        if isempty(curr_data)
            continue;
        end
        
        timestamps(j) = mean(wifi_timestamps(unique_data(:,1)==i-1));
        [~,inds] = sort(abs(img_timestamps - timestamps(j)));
        inds = inds(1:2);
        times = img_timestamps(inds);
        if (times(2)-timestamps(j))*(times(2)-timestamps(j)) > 0
            c1_p_ci = g_ba_ci(:,inds(1)-init_ba_indx+1);
        else
            c1_p_ci = (g_ba_ci(:,inds(1)-init_ba_indx+1) * abs(times(2)-timestamps(j)) + g_ba_ci(:,inds(2)-init_ba_indx+1) * abs(times(1)-timestamps(j))) / abs(times(2)-times(1));
        end
        pt_xy = c1_p_ci([1,2]);     %xy
        b(2*j-1:2*j) = pt_xy;
        
        sig = zeros(1, Ns);
        sig(curr_data(:,1)) = curr_data(:,2);
        sig_normalized = sig / norm(sig);
        
        A(2*j-1, 1:Ns) = sig_normalized;
        A(2*j, Ns+1:end) = sig_normalized;
        j = j+1;
    end
    
    N = j-1;
    A(2*N+1:end, :) = [];
    b(2*N+1:end) = [];
    timestamps(N+1:end) = [];

    start_ind = test_data_num*2+1;
    x = (A(start_ind:end,:)'*A(start_ind:end,:))\(A(start_ind:end,:)'*b(start_ind:end));
    % model = reshape(x, [2 length(x)/2]);
    figure, axis equal, 
    plot(g_ba_ci(1,:), g_ba_ci(2,:),'r')
    
    for i=1:N
        pt = A(2*i-1:2*i, :)*x;
        hold on, plot(pt(1),pt(2),'go')
        hold on, plot(b(2*i-1),b(2*i),'bx')
        plot([pt(1), b(2*i-1)], [pt(2), b(2*i)])
    end
    
    sig_mat = A(1:2:end, 1:Ns);
    sig_pos = reshape(b,2,N);
    
%     cell_size = 1;
%     pts = g_ba_ci(1:2, :);
%     pt_x = b(1:2:end);
%     pt_y = b(2:2:end);
%     min_x = min([pts(1,:), pt_x']);
%     min_y = min([pts(2,:), pt_y']);
%     
%     pts(1,:) = pts(1,:)-min_x;
%     pts(2,:) = pts(2,:)-min_y;
%     pt_x = pt_x-min_x;
%     pt_y = pt_y-min_y;
%     temp_grid = zeros(ceil(cell_size*max([pts(1,:), pt_x']))+1, ceil(cell_size*max([pts(2,:), pt_y']))+1);
%     
%     st_ind = 38;
%     st_levels = A(1:2:end,st_ind);
% 
%     for i=1:length(pt_x)
%         temp_grid(ceil(cell_size*pt_x(i))+1, ceil(cell_size*pt_y(i))+1) = st_levels(i);
%     end
%     
%     figure,mesh(temp_grid)
    
else
    A = zeros(3*N, 3*Ns);
    b = zeros(3*N, 1);
    timestamps = zeros(N, 1);
    j = 1;
    
    for i = 1:N
        % curr_data = unique_data(unique_data(:,1)==i-1, [2,3]);curr_data(:,2)=-curr_data(:,2);
        curr_data = unique_data(unique_data(:,1)==i-1, [2,4]);
        if isempty(curr_data)
            continue;
        end
        
        timestamps(j) = mean(wifi_timestamps(unique_data(:,1)==i-1));
        [~,inds] = sort(abs(img_timestamps - timestamps(j)));
        inds = inds(1:2);
        times = img_timestamps(inds);
        c1_p_ci = (g_ba_ci(:,inds(1)-init_ba_indx+1) * abs(times(2)-timestamps(j)) + g_ba_ci(:,inds(2)-init_ba_indx+1) * abs(times(1)-timestamps(j))) / abs(times(2)-times(1));
        pt_xy = c1_p_ci([1,2, 3]);     %xyz
        b(3*j-2:3*j) = pt_xy;
        
        sig = zeros(1, Ns);
        sig(curr_data(:,1)) = curr_data(:,2);
        sig_normalized = sig / norm(sig);
        
        A(3*j-2, 1:Ns) = sig_normalized;
        A(3*j-1, Ns+1:2*Ns) = sig_normalized;
        A(3*j, 2*Ns+1:3*Ns) = sig_normalized;
        j = j+1;
    end
    
    N = j-1;
    A(3*N+1:end, :) = [];
    b(3*N+1:end) = [];
    timestamps(N+1:end) = [];

    start_ind = test_data_num*3+1;
    x = (A(start_ind:end,:)'*A(start_ind:end,:))\(A(start_ind:end,:)'*b(start_ind:end));
    % model = reshape(x, [3 length(x)/3]);
    plot3(g_ba_ci(1,:), g_ba_ci(2,:), g_ba_ci(3,:),'r')
    
    for i=1:test_data_num
        pt = A(3*i-2:3*i, :)*x;
        hold on, plot3(pt(1),pt(2), pt(3),'go')
        plot3(b(3*i-2),b(3*i-1),b(3*i),'bx')
    end
end
