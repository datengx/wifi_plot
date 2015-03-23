clear all;
% close all;

% if exist('nameAndIdx.mat')
%     delete('nameAndIdx.mat');
% end
% if exist('wifiStationAndIdx.mat')
%     delete('wifiStationAndIdx.mat');
% end

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

% unique_data(unique_data(:,3)<=-85,:) = [];
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
            if freqs(l)>5000
                continue;
            end
            freq_data = ap_data(ap_data(:,6)==freqs(l), :);
            mean_freq_data = mean(freq_data,1);
            temp_data = [temp_data; mean_freq_data];
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
idAndPosition = [5,  11.48,  -4.487,  1.2192;...
                 6, -20.531,  10.533, 1.2192;...
                10, -13.887, -16.508, 1.1684;...
                14, -32.955, -16.996, 1.27;...
                20, -3.84,   -24.83,  1.1684;...
                21, -12.17,  -47.06,  1.4732;...
                22, -28.98,  -45.769, 1.4982];
            
%% Least square to solve n and A
[N,M] = size(unique_data);
unique_data = [unique_data,zeros(N,3)];
% z_s_avg = unique_data(:,8)/N;
for k = 1:N
    unique_data(k,9:11) = idAndPosition(...
                         idAndPosition(:,1) == unique_data(k,2),[2,3,4]);
end

R = unique_data(:,3);
D = 10*log10( sqrt((unique_data(:,6)-unique_data(:,9)).^2 ...
                 + (unique_data(:,7)-unique_data(:,10)).^2 ...
                 + (unique_data(:,8)-unique_data(:,11)).^2)...
               );

n_A = [D,-ones(length(R),1)] \ (-R-20*log10(unique_data(:,4)))
% [D,ind] = sort(D);
% plot(D,-n_A(1)*D+n_A(2)-20*log10(mean(unique_data(ind,4))))
% hold on,plot(D,R(ind),'.')


ap_ids = unique(unique_data(:,2));
colors = jet(length(ap_ids));
for i = 1%:length(basement_ap_ids)
    figure, axis equal, hold on
    plot3(g_ba_ci(1,:), g_ba_ci(2,:), g_ba_ci(3,:),'r')
    %% Plot WIFI position on trajectory
    plot3(idAndPosition(:,2),idAndPosition(:,3),idAndPosition(:,4),'go')

    j = ap_ids(i);
    curr_data = unique_data(unique_data(:,2)==j,:);
    freqs = unique(curr_data(:,4));
    for k = 1:size(curr_data,1)
        d = 10.^((n_A(2)-curr_data(k,3)-20*log10(curr_data(k,4)))/10/n_A(1));
        circle(curr_data(k,6), curr_data(k,7), sqrt(d^2 - curr_data(k,8)^2));
    end
    scatter3(curr_data(:,6), curr_data(:,7), curr_data(:,8)+2*i, 'MarkerEdgeColor', colors(i,:), 'Marker', 'x')

    figure, hold on, axis equal
    for l = 1:length(freqs)
        scatter3(curr_data(curr_data(:,4)==freqs(l),6),curr_data(curr_data(:,4)==freqs(l),7),curr_data(curr_data(:,4)==freqs(l),3),'MarkerEdgeColor', colors(l*floor(size(colors,1)/2),:), 'Marker', '.')
    end
end
% 
% % figure, axis equal
% for i = 1:length(basement_ap_ids)
%     j = ap_ids(i);
%     disp(['station', num2str(j)])
%     curr_data = unique_data(unique_data(:,2)==j,:);
%     rssi = R(unique_data(:,2)==j,:);
%     n_A_pre = n_A;
% %     for k=1:2
%         d = 10.^((n_A_pre(2)-rssi-20*log10(curr_data(:,4))) / 10 / n_A_pre(1));
%         
%         a=curr_data(:,6).*curr_data(:,6); a = a - mean(a);
%         b=curr_data(:,7).*curr_data(:,7); b = b - mean(b);
%         c = d.*d; c = c - mean(c);
%         e = (a+b-c)/2;
%         
%         A = [curr_data(:,6)-mean(curr_data(:,6)), curr_data(:,7)-mean(curr_data(:,7))];
%         x = A\e;
%         
%         z = d.*d - (curr_data(:,6)-x(1)).^2 - (curr_data(:,7)-x(2)).^2;
%         disp([sum(z<=0), sum(z>0)])
%         z = z(z>0);
%         z = sqrt(mean(z)) - mean(curr_data(:,8));
%         disp([x', z])
%         
% %         D = 10*log10( sqrt((curr_data(:,6)-x(1)).^2 ...
% %             + (curr_data(:,7)-x(2)).^2 ...
% %             + (curr_data(:,8)-curr_data(:,11)).^2)...
% %             );
% %         
% %         n_A_pre = [D,-ones(length(D),1)] \ (-rssi-20*log10(curr_data(:,4)))
% % %         n_A_pre(2) = mean(rssi+n_A_pre(1)*D+20*log10(curr_data(:,4)));
% % %         [D,ind] = sort(D);
% % %         plot(D,-n_A_pre(1)*D+n_A_pre(2)-20*log10(mean(curr_data(:,4))))
% % %         hold on,plot(D,rssi(ind),'.')
% %     end
% end

% figure, axis equal, hold on
% plot(g_ba_ci(1,:), g_ba_ci(2,:),'r')
% plot(idAndPosition(:,2),idAndPosition(:,3),'go')
XYZ=[];
for i = 1:length(basement_ap_ids)
    j = ap_ids(i);
    disp(['station', num2str(j)])
    curr_data = unique_data(unique_data(:,2)==j,:);
    rssi = unique_data(unique_data(:,2)==j,3);
    n_A_pre = n_A;

    d = 10.^((n_A_pre(2)-rssi-20*log10(curr_data(:,4))) / 10 / n_A_pre(1));
    
    a=curr_data(:,6).*curr_data(:,6); a = a - mean(a);
    b=curr_data(:,7).*curr_data(:,7); b = b - mean(b);
    c = d.*d; c = c - mean(c);
    e = (a+b-c)/2;
    A = [curr_data(:,6)-mean(curr_data(:,6)), curr_data(:,7)-mean(curr_data(:,7))];
    x = A\e;
    
    z = 1.2;
    xyz = [x;z];
    err = 1e9;
    prev_err = 1e12;
%     del = .3;
    while err>1e-3 && (prev_err-err)>1e-3
        J = [2*(x(1)-curr_data(:,6)), 2*(x(2)-curr_data(:,7)), 2*xyz(3)*ones(size(curr_data,1),1)];
        r = (x(1)-curr_data(:,6)).^2 + (x(2)-curr_data(:,7)).^2 + xyz(3)^2 - (d).^2;
%         w = del*(sqrt(1+(r/del).^2)-1);
%         W=diag(w);
%         dxyz = (J'*W*J)\(J'*W*r);
        dxyz = (J'*J)\(J'*r);
        xyz_prev = xyz;
        xyz = xyz+dxyz;
        prev_err = err;
        err = mean(r.*r);
    end
    if prev_err<err
        xyz=xyz_prev;
    end
    disp(xyz')
    XYZ = [XYZ, xyz];   

end

% plot(XYZ(1,:),XYZ(2,:),'bx')
% title('AP locations')
% xlabel('x'), ylabel('y')
% legend('Trajectory','Ground truth AP locations','Estimated AP locations')