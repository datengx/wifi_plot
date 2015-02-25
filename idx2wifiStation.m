function [ wifiStation ] = idx2wifiStation( idx )
%idx2wifiStation 
% 

    load('wifiStationAndIdx.mat');
    assert(size(idx,2)==1);
    wifiStation = [];
    for j = 1:size(idx,1)
        i = find(wifiStationAndIdx(:,2)==idx(j,1));
        if size(i, 1) ~= 0
            wifiStation = [wifiStation; wifiStationAndIdx(i,1)];
        else
            wifiStation = [wifiStation; -1];
        end
    end
end

