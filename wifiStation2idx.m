function [ idx ] = wifiStation2idx( wifiStation )
%wifiStation2idx convert from wifiStation2idx
%
%   Convert from wifistation name to a unique idx
%
assert(size(wifiStation,2)==1);
if exist('wifiStationAndIdx.mat')
    load('wifiStationAndIdx.mat');
else
    wifiStationAndIdx= [];
end
idx = [];

for j = 1:size(wifiStation,1)
    wifiStation(j,1) = bitand(uint64(wifiStation(j,1)), uint64(hex2dec('fffffffffff0')), 'uint64');
    if ~isempty(wifiStationAndIdx)
%         i = find(wifiStationAndIdx(:,1)==wifiStation(j,1));
        i = find(abs(wifiStationAndIdx(:,1)-wifiStation(j,1))<=16);
    else
        i = [];
    end
    if size(i, 1) ~= 0
        idx = [idx; wifiStationAndIdx(i,2)];
    else
        wifiStationAndIdx(size(wifiStationAndIdx, 1)+1,:) = [wifiStation(j,1), size(wifiStationAndIdx, 1)+1];
        idx = [idx; size(wifiStationAndIdx, 1)];
    end
end
save('wifiStationAndIdx.mat', 'wifiStationAndIdx');

end

