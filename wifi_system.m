clear all;
close all;

path = '/home/da/dev/wifi_plot/mrinal_wifi_dataset/wifi/';
files = dir(path);
%% The first two are '.' and '..'
file_names = {files(3:end).name}';

data = [];
id_names = {};
for i = 1:length(file_names)
    fid = fopen([path, file_names{i}]);
    token = textscan(fid, '%d64%d64%d64%f64%d64%s\n');
       
    data = [data; double((i-1)*ones(length(token{1}), 1)), wifiStation2idx(double(token{1})), double(token{2}), (10.^(double(token{2})./10))*1000, double(token{3}), double(token{4}), double(token{5}), double(name2idx(token{6}))];
%     id_names{length(id_names):length(data)} = token{6};
    fclose(fid);
end

%% remove duplicate
[C, IA, IC] = unique(data(:,[2, 5]), 'rows', 'stable');
%% create data matrix without duplicate
data = data(IA, :);


%%
%   WIFI data processing: Plot all signal level of all the location per
%   station
%   Compare WIFI data collecting from linux and yellowstone
%   MARS
%   06/02/2015
%
%%

clear all;

%% add dataset path
DATAPATH = ['./spot01'; './spot02'; './spot03'; './spot04'; './spot05'; './spot06'; './spot07'; './spot08'; './spot09'; './spot10'];
file_size = [];
file_size_ys = [];
for text_idx = 1:size(DATAPATH, 1)
    timestamp = dlmread(strcat(DATAPATH(text_idx,:), '/linux/timestamp.txt'));
    file_size = [file_size, size(timestamp, 1)];
    cd(strcat(DATAPATH(text_idx, :), '/ys'));
    d = dir;
    cd('../..');
    file_size_ys = [file_size_ys, size(d, 1) - 2];
end


data = [];
data_ys = [];

%% Extract linux data
%% Create a dataset
dataset = cell(size(DATAPATH, 1),1);
dataset_ys = cell(size(DATAPATH, 1),1);
for text_idx = 1:size(DATAPATH, 1)
    data = [];
    for i = 1:file_size(text_idx)-1
       filename = sprintf('%s/linux/%08d.txt', DATAPATH(text_idx, :), i);
       token = dlmread(filename);
       data = [data; token];
    end
    dataset{text_idx, 1} = data;
end
%% Extract yellowstone data
for text_idx = 1:size(DATAPATH, 1)
    data_ys = [];
    for i = 1:file_size_ys(text_idx)-1
       filename_ys = sprintf('%s/ys/%08d.txt', DATAPATH(text_idx, :), i);
       fid_ys = fopen(filename_ys);
       token_ys = textscan(fid_ys, '%d64%d64%d64%f64%d64%s\n');
       data_ys = [data_ys; token_ys{1}, token_ys{2}, token_ys{3}, token_ys{4}, token_ys{5}];
       fclose(fid_ys);
    end
    dataset_ys{text_idx, 1} = data_ys;
end
%% Get linux station
station = cell(size(DATAPATH, 1),1);
station_ys = cell(size(DATAPATH, 1),1);
for text_idx = 1:size(DATAPATH, 1)
    data = dataset{text_idx, 1};
    data = data(:,1);
    data = unique(data);
    station{text_idx, 1} = data(2:end);
end
%% Get yellowstone station
for text_idx = 1:size(DATAPATH, 1)
    data_ys = dataset_ys{text_idx, 1};
    data_ys = data_ys(:,1);
    data_ys = unique(data_ys); % unique station
    station_ys{text_idx, 1} = data_ys(1:end);
end
%% Find out the common station
common_station = station{1, 1};
common_station_ys = station_ys{1, 1};
for text_idx = 1:size(DATAPATH, 1)-1
    common_station = union(common_station, station{text_idx+1, 1});
    common_station_ys = union(common_station_ys, station_ys{text_idx+1, 1});
end

common_station_overall = union(common_station, common_station_ys);
%% label all the station
len = length(common_station_overall);
station_label = (1:len)';
station_overall_and_label = [common_station_overall, station_label]; % station (first column) with label (second column)
save('station_overall_and_label.mat', 'station_overall_and_label');
load('station_overall_and_label.mat');

full_data = [];
for x = 1:size(station_overall_and_label, 1)
    overall_data = [];
    overall_data_ys = [];
    station_find = common_station_overall(x);
    station_find_ys = station_find;

    %% Calculate average signal level
    for text_idx = 1:size(DATAPATH, 1)
        %% Creating plot for linux
        data = dataset{text_idx, 1};
        idx = find(data(:,1) == station_find);
        data = data(idx,:);
        data = data(:,2);
        data = data(data<0);
        data = data(data>-100);

        %% Creating plot for yellowstone
        data_ys = dataset_ys{text_idx, 1};
        idx_ys = find(data_ys(:,1) == station_find_ys);
        data_ys = data_ys(idx_ys, :);
        data_ys = data_ys(:,2);
        data_ys = data_ys(data_ys<0);
        data_ys = data_ys(data_ys>-100);

        %% linux WIFI level signal
        overall_data = [overall_data; [ones(length(data),1)*text_idx, data]];
        %% yellowstone level signal
        overall_data_ys = [overall_data_ys; [ones(length(data_ys),1)*text_idx, data_ys]];

    end
    
    station_vec = ones(length(overall_data), 1)*double(returnLabelWithStation(station_find));
    station_vec_ys = ones(length(overall_data_ys), 1)*double(returnLabelWithStation(station_find_ys));
%     full_data = [full_data; station_vec, overall_data];
    full_data = [full_data; station_vec_ys, overall_data_ys];
end
  
st_1st = min(full_data(:,1));
st_last = max(full_data(:,1));
sig_normalized = zeros(max(full_data(:,1)), max(full_data(:,2))-min(full_data(:,2))+1);
st_len = zeros(10,1);
for i=min(full_data(:,2)):max(full_data(:,2))
    sig = zeros(st_last, 1);
    ind = full_data(:,2)==i;
    dd = full_data(ind, [1,3]);
%     [~,index] = sort(dd(:,1));
%     dd= dd(index,:);
    stations = unique(dd(:,1));
    st_len(i)=length(stations);
    for j = 1:length(stations)
        strengths = dd(dd(:,1)==stations(j), 2);
        avg_strength = mean(strengths);
        sig(stations(j)) = avg_strength;
    end
    sig_normalized(:,i) = sig / norm(sig);
end

ss = zeros(max(full_data(:,2))-min(full_data(:,2))+1);
for i=min(full_data(:,2)):max(full_data(:,2))
    for j=min(full_data(:,2)):max(full_data(:,2))
        ss(i, j) = sig_normalized(:,i)' * sig_normalized(:,j);
    end
end

ssp = ss - eye(max(full_data(:,2))-min(full_data(:,2))+1);
rank(ss)
eigs = flip(eig(ss));
plot(eigs)

