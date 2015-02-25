function [ idx ] = name2idx( name )
%name2idx 
%   
    load('nameAndIdx.mat');
    idx = [];
    for j = 1:size(name,1)
        i = 0;
        sz_current = size(nameAndIdx, 1);
        for k = 1:sz_current
            i = strcmp(nameAndIdx{k,1},name{j,1});
            if i ~= 0
                idx = [idx; nameAndIdx{k,2}];
                break;
            end
        end
        
        if i == 0
            nameAndIdx{sz_current+1,1} = name{j,1};
            nameAndIdx{sz_current+1,2} = sz_current+1;
            idx = [idx; sz_current+1];
        end
    end
    save('nameAndIdx.mat', 'nameAndIdx');
end

