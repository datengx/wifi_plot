function [ idx ] = wifiStation2idx( wifiStation )
%wifiStation2idx convert from wifiStation2idx
%   
%   Convert from wifistation name to a unique idx
%
    assert(size(wifiStation,2)==1);
    load('wifiStationAndIdx.mat');
    idx = [];
    for j = 1:size(wifiStation,1)
        i = find(wifiStationAndIdx(:,1)==wifiStation(j,1));
        if size(i, 1) ~= 0
            idx = [idx; wifiStationAndIdx(i,2)];
        else
            wifiStationAndIdx(size(wifiStationAndIdx, 1)+1,:) = [wifiStation(j,1), size(wifiStationAndIdx, 1)+1];
            idx = [idx; size(wifiStationAndIdx, 1)];
        end
    end
    save('wifiStationAndIdx.mat', 'wifiStationAndIdx');
end
