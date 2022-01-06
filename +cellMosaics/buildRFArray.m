% Builds a spatial array of cells using a defined receptive field and
% cell spacing (in microns)

function [mosaic,centerCoordinates,cellTracker] = buildRFArray(obj, RFFilter)
    
    if strcmp(obj.rfRegion,'center')
        
        % Crop RF to center
        [x,y] = find(RFFilter > 0);
        RF_cropped = RFFilter(min(x):max(x),min(y):max(y));
        RF_cropped = RF_cropped .* (RF_cropped > 0); % Remove all surround

    else
        
        %%% Crop RF to space that captures 95% of inhibitory behavior
        [xx,yy] = meshgrid(1:obj.videoSize(2),1:obj.videoSize(1));
        r = sqrt((xx - obj.videoSize(2)/2).^2 + (yy - obj.videoSize(1)/2).^2); 
        
        % Calculate cumulative inhibiton for different radii
        total_inhibition = sum(RFFilter(RFFilter < 0));
        incremental_inhibition = zeros(max(size(RFFilter)),1);
        for a = 1:max(size(RFFilter))
            tmp = RFFilter .* (r <= a);
            incremental_inhibition(a) = sum(tmp(tmp < 0));
        end
        
        % Crop to 95% inhibition
        [~,i] = min(abs(incremental_inhibition - total_inhibition*0.95));
        [x,y] = find(r < i);
        RF_cropped = RFFilter(min(x):max(x),min(y):max(y));
        
        if strcmp(obj.rfRegion,'surround')
            RF_cropped = RF_cropped .* (RF_cropped < 0); % Only include surround
        elseif strcmp(obj.rfRegion,'full-field')
            RF_cropped = cat(3,RF_cropped .* (RF_cropped > 0), RF_cropped .* (RF_cropped < 0));
        end
    end

    % To view cropped receptive-field:
%     mesh(sum(RF_cropped,3))
    
    % Divide into wedge-shaped radial slices (as needed)
    if obj.slices > 1
        [xx,yy] = meshgrid(1:size(RF_cropped,2),1:size(RF_cropped,1));
        
        % Radial space
        th = atan((xx - size(RF_cropped,2)/2) ./ (yy - size(RF_cropped,1)/2));
        th = abs(th-pi/2);
        nonsmooth = find(diff(th) > pi/2,1);
        th(1:nonsmooth,:) = th(1:nonsmooth,:) + pi;
        th = rad2deg(th);
        
        % Divide into slices
        theta = 0:360/obj.slices:360;
        RF_cropped_sliced = zeros([size(RF_cropped,1), size(RF_cropped,2), obj.slices.*size(RF_cropped,3)]);
        for c = 1:length(theta)-1
            angFilt = (th >= theta(c)) & (th < theta(c+1)); % Angular filter (theta)
            RF_cropped_sliced(:,:,c) = RF_cropped(:,:,1) .* angFilt;
            
            % For full-field stimuli, include the surround
            if size(RF_cropped,3) > 1
                RF_cropped_sliced(:,:,c+obj.slices) = RF_cropped(:,:,2) .* angFilt;
            end
        end
    else
        RF_cropped_sliced = RF_cropped;
    end
    
    % This produces an n x m x s array, where the third dimension (s)
    % represents each isolated region of the receptive-field
    
    % To view cropped and sliced receptive-field:
%     for a = 1:size(RF_cropped_sliced,3)
%         imagesc(RF_cropped_sliced(:,:,a))
%         waitforbuttonpress % Press any button to view next segment
%     end

    % Convert all variables to arcmin
    sz = size(RF_cropped_sliced,1); % Diameter of individual RGCs (already in arcmin)
    sp = round(cellMosaics.utils.changeUnits(obj.spacing,obj.micronsPerPixel,'um2arcmin')); % Spacing
    n = round(cellMosaics.utils.changeUnits(obj.noise,obj.micronsPerPixel,'um2arcmin')); % Noise
    n = n ./ sqrt(2); % Maximum magnitude = expected noise

    % As we tile our cells, RGCs that lie on the edge of our "retina"
    % encode differently than cells that lie in the center.
    
    % To combat edge effects, we will create a larger retina.
    % We will tile cells in our larger retina and then crop to size.
    
    % Create array of cells
    videoSize_larger = round(obj.videoSize + sp*10); % Identify larger space
    
    % Identify coordinates
    x = floor(sz/2):sp:(videoSize_larger(1)-ceil(sz/2));
    y = floor(sz/2):sp:(videoSize_larger(2)-ceil(sz/2));
    if strcmp(obj.tilingStyle,'hex')
        x = floor(sz/2):round(sp .* (sqrt(3)/2)):(videoSize_larger(1)-ceil(sz/2));
    end

    mosaic = zeros([videoSize_larger, length(y) .* length(x) .* size(RF_cropped_sliced,3)]);
    centerCoordinates = zeros(length(y) .* length(x) .* size(RF_cropped_sliced,3),2);
    cellTracker = zeros(length(y) .* length(x) .* size(RF_cropped_sliced,3),1);
    counter = 1;
    for a = 1:length(x)
        for b = 1:length(y)
            
            % Define coordinates
            x_coordinate = x(a) + round(rand()*n);
            y_coordinate = y(b) + round(rand()*n);
            if mod(a,2) == 0 & strcmp(obj.tilingStyle,'hex')
                y_coordinate = y_coordinate + round(sp / 2); % Offset for hexagonal tiling
            end
                            
            % Map out full range of coordinates for placing receptive field
            xx = round(x_coordinate - sz/2) : (round(x_coordinate - sz/2) + sz - 1);
            yy = round(y_coordinate - sz/2) : (round(y_coordinate - sz/2) + sz - 1);

            % Check coordinates are valid
            if sum([xx,yy] <= 0) == 0 && sum(xx > videoSize_larger(1)) == 0 && ...
               sum(yy > videoSize_larger(2)) == 0      
                
                % Place receptive-field
                for c = 1:size(RF_cropped_sliced,3)
                    tmp = zeros(videoSize_larger);
                    tmp(xx,yy) = tmp(xx,yy) + RF_cropped_sliced(:,:,c);
                    
                    % Retina mosaic for export
                    mosaic(:,:,counter) = mosaic(:,:,counter) + tmp;
                    
                    % Save coordinates to check distance
                    centerCoordinates(counter,:) = [x_coordinate, y_coordinate];
                    
                    % Keep track of which sliced receptive-fields belong to
                    % each cell
                    cellTracker(counter) = counter - c + 1;
                    counter = counter+1;
                end
            end
        end
    end
    
    % Crop s.t. top left neuron lies at the top left corner
    [xx,yy] = find(sum(mosaic,3) > 0);
    offset = [min(xx) min(yy)] + round(sp*2);
    mosaic = mosaic(offset(1):offset(1)+obj.videoSize(1)-1,offset(2):offset(2)+obj.videoSize(2)-1,:);
    centerCoordinates = cellMosaics.utils.changeUnits(centerCoordinates,obj.micronsPerPixel,'arcmin2um'); % Spacing
    [~,~,cellTracker] = unique(cellTracker);
    
    % Remove unplaced neurons (that likely impinge on edges of space)
    mosaic = mosaic(:,:,centerCoordinates(:,1) > 0);
    cellTracker = cellTracker(centerCoordinates(:,1) > 0);
    centerCoordinates = centerCoordinates(centerCoordinates(:,1) > 0,:);
    [~,~,cellTracker] = unique(cellTracker,'rows');
end

