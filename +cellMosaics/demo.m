%%%
% This script builds a mosaic of retinal ganglion cell (RGC) receptive
% fields. The resulting collection of RGCs can be used to simulate how 
% populations of RGCs might spatially integrate stimuli in the retina.
% by J. Freedland, 12/2021
%%%

% Select cell type to place
clear
cellType = 'on parasol';
% cellType = 'on midget';

%% Generate receptive fields for single neurons
if strcmp(cellType,'on parasol')
    
    % Classical receptive fields use a difference-of-gaussians fit.
        % centerSigma: sigma of gaussian in receptive-field center.
        % surroundSigma: sigma of gaussian in receptive-field surround.
    obj.centerSigma   = 50; % in microns
    obj.surroundSigma = 160; % in microns
    
    % Define spacing of neurons in the retina
    obj.spacing = 150; % in microns
    
elseif strcmp(cellType,'on midget')

    obj.centerSigma   = 30; % in microns
    obj.surroundSigma = 100; % in microns
    obj.spacing = 75; % in microns

end

% Create receptive field
% Output units: 1 pixel = 1 arcmin
obj.monitorSize         = [600 800]; % Size of monitor [height, width]
obj.micronsPerPixel     = 1.65;      % How many microns each pixel spans
obj.videoSize = cellMosaics.utils.changeUnits(obj.monitorSize,obj.micronsPerPixel,'pix2arcmin');
RFFilter = cellMosaics.utils.calculateFilter(obj);
mesh(RFFilter)

%% We now tile receptive fields into an array
% THIS SECTION MAY TAKE UP TO 1 MINUTE TO RUN - speed increases in development

% Portion of each neuron's receptive-field to tile:
%   'center': only places the receptive-field center of each neuron
%   'surround': only places the receptive-field surround of each neuron
%   'full-field': places the entire receptive-field
obj.rfRegion = 'full-field';

% Tiling style: How to arrange neurons in our grid:
%   'hex': hexagonal closest packed
%   'grid': simple grid
obj.tilingStyle = 'hex';

% Neurons aren't perfectly X microns apart.
    % Noise offsets each cell by a random value between [0 obj.noise]
    % in x, y directions
obj.noise = obj.spacing/4; % in microns

% Based on a paper by J. Freedland & F. Rieke (2021), dividing a receptive field 
% into 8 "pie slices" significantly improves prediction. Set to 1 to
% place a uniform center/surround.
obj.slices = 8;

% Build mosaic of neurons
% Models with no slices (set obj.slices = 1) run much quicker
tic
[mosaic, coordinates, cellTracker] = cellMosaics.buildRFArray(obj, RFFilter);
toc
imagesc(sum(mosaic(:,:,1:end),3))

%% Check cells are proper distance apart
averageDistance = round(cellMosaics.utils.checkCellSpacing(coordinates)); % in microns
disp(strcat('Actual distance between RGCs (on average): ',mat2str(averageDistance),' um'))

%% Create a stimulus
% Grating width
width = 200; % in microns
for w = 1:length(width)
    
    gratingWidth = round(cellMosaics.utils.changeUnits(width(w),obj.micronsPerPixel,'um2arcmin')); % Convert to arcmin

    % Create grating
    % Note: stimulus should ALWAYS be in units of contrast (typically -1 to 1)
    stimulus = repmat(repelem([0.9 -0.9],1,gratingWidth),size(mosaic,1),size(mosaic,2));
    stimulus = stimulus(1:size(mosaic,1),1:size(mosaic,2));
end
imagesc(stimulus)

%% Spatially integrate stimulus using mosaic
% Linearly integrate stimulus into each receptive-field section
% Note: all NaN values lie outside of mosaic array
linearValues = cellMosaics.utils.linearIntegration(mosaic,cellTracker,stimulus);

% Input values into nonlinear model
if strcmp(obj.rfRegion, 'center')
    centerValues = linearValues(:,1:obj.slices);
    surroundValues = zeros(size(centerValues));
elseif strcmp(obj.rfRegion, 'surround')
    surroundValues = linearValues(:,1:obj.slices);
    centerValues = zeros(size(surroundValues));
elseif strcmp(obj.rfRegion, 'full-field')
    centerValues = linearValues(:,1:obj.slices);
    surroundValues = linearValues(obj.slices+1:end);
end

% Note: all NaN values lie outside of mosaic array
responses = cellMosaics.utils.spatialModelFitted(centerValues,surroundValues,cellType);
realResponses = responses(~isnan(responses)); % Ignores edge effects

%% Visualize all aspects together
responseMosaic = cellMosaics.utils.buildResponseMosaic(mosaic,cellTracker,responses);

figure(1)
subplot(1,3,1)
imagesc(stimulus)
title('stimulus')
subplot(1,3,2)
imagesc(sum(mosaic,3))
title('cell mosaic')
subplot(1,3,3)
imagesc(responseMosaic)
title('cell responses')