% Uses each cell in a mosaic to calculate a stimulus' linearly integrated
% value
function linearValues = linearIntegration(mosaic,cellTracker,stimulus)
    
    % Stimulus should have input units as contrast (typically, -1 to 1)
    
    n = max(cellTracker);
    m = sum(cellTracker == 1);
    linearValues = zeros(n,m);

    for a = 1:max(unique(cellTracker))

        % Identify individual cell
        M = mosaic(:,:,cellTracker == a);

        % Calculate linear equivalent disk
        for b = 1:size(M,3)
            A = sum(M(:,:,b) .* stimulus, [1 2]); % RF .* stimulus .* region
            B = sum(M(:,:,b), [1 2]); % RF .* region
            linearValues(a,b) = (A / B);
        end
    end
end

