function responseMosaic = buildResponseMosaic(mosaic,cellTracker,responses)

    responseMosaic = zeros(size(mosaic,1),size(mosaic,2));
    
    for a = 1:max(unique(cellTracker))
        
        % Identify individual cell
        if ~isnan(responses(a))
            M = sum(mosaic(:,:,cellTracker == a),3);

            % Place response value in cell area
            responseMosaic = responseMosaic + (M .* responses(a));
        end
    end
end

