%%%
% Converts spatial luminances into spike outputs.
% JMF 06/2021
%
% INPUTS: centerLuminances: N x M matrix of M regional luminances across N images
%         surroundLuminances: N x M matrix of M regional luminances across N images
%
%              Luminances should be inputed as fractional contrasts,
%              defined as: (luminance - backgroundLuminance) / backgroundLuminance
%
%              NOTE: The ordering of luminances (location of each region) in centerLuminances and 
%              surroundLuminance have no effect on the output of the model.
%
%         type: whether the cell is an 'on' or 'off' pathway cell. This
%         model uses empirically measured parameters from both pathways.
%
% OUTPUTS: sp: relative number of spikes predicted for each image.
%
%               i.e sp = 0.3 implies that flashing an image
%               will result in 30% of the cell's maximum spiking rate.
%%%

function sp = spatialModelFitted(centerLuminances,surroundLuminances,type)

    centerInput     = centerLuminances;
    surroundInput   = surroundLuminances;
    
    % If different numbers of center and surround regions,
    % create degenerate dimensions.
    if size(centerInput,2) ~= size(surroundInput,2)
        reps = lcm(size(centerInput,2),size(surroundInput,2));
        centerInput = repelem(centerInput,1,reps/size(centerInput,2));
        surroundInput = repelem(surroundInput,1,reps/size(surroundInput,2));
    end

    % Empirically-defined parameters
    if contains(type,'on parasol')
        params = [0.772561921672811,-0.544239686227808,-0.146423253733880,0.854487003405620];
    elseif contains(type,'off parasol')
        params = [0.699063417620684,-0.418676713678296,-0.215891795996638,1.04047440129262];
    elseif contains(type, 'on midget')
        params = [0.619208019, -0.618048319, -0.040699022, 4.678914178];
    elseif contains(type,'off midget')
        params = [0.585648295, -0.313180092, -0.010224563, 3.890761152];
    end
    
    % Indices for parameters
    centerLinearIndex       = 1;
    surroundLinearIndex     = 2;
    centerNonlinearIndex    = 3:4;
    
    %%% Spatial model
    % This code originally used a gradient descent to fit parameters, and has been
    % retrofitted to allow for static parameters.
    model.center       = @(params) centerInput; % Raw center
    model.surround     = @(params) surroundInput; % Raw surround

    % Weight center and surround uniformly, then combine.
    model.combined     = @(params) model.center(params) .* params(centerLinearIndex) + model.surround(params) .* params(surroundLinearIndex);
    
    % Rectify
    model.negativeRectifier    = @(params) (model.combined(params) > params(centerNonlinearIndex(1))) .* model.combined(params) + ...
                                           (model.combined(params) <= params(centerNonlinearIndex(1))) .* params(centerNonlinearIndex(1));
    model.positiveRectifier    = @(params) (model.negativeRectifier(params) < params(centerNonlinearIndex(2))) .* model.negativeRectifier(params) +...
                                           (model.negativeRectifier(params) >= params(centerNonlinearIndex(2))) .* params(centerNonlinearIndex(2));
    
    % Integrate
    model.linearity            = @(params) nanmean(model.positiveRectifier(params),2);
    sp = model.linearity(params);
end

