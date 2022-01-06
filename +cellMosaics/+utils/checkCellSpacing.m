% Checks cell spacing after tiling mosaic

function averageDistance = checkCellSpacing(coordinates)
    
    uniqueCoords = unique(coordinates,'rows');
    averageDistance = zeros(size(uniqueCoords,1),1);
    
    for a = 1:size(uniqueCoords,1)
        
        % Calculate euclidean distance
        d = sqrt(sum((uniqueCoords(a,:) - uniqueCoords).^2,2)); % Euclidean distance
        d = sort(d,'ascend'); % in microns
        d = d(2:end); % first value is always 0 (cell's distance from itself)
        
        % Check four closest nearest neighbors
        averageDistance(a) = nanmean(d(1:4));
    end
    
    averageDistance = nanmean(averageDistance); % Average across all cells
end

