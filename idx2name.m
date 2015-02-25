function [ name ] = idx2name( idx )
%idx2name 
%   Detailed explanation goes here
    load('nameAndIdx.mat');
    name = cell(0,1);
    sz = size(nameAndIdx,1);
    for i = 1:size(idx,1)
        for j = 1:sz
            if idx(i,1) == nameAndIdx{j,2}
                name{size(name, 1)+1, 1} = nameAndIdx{j,1};
            end
        end
    end
end

