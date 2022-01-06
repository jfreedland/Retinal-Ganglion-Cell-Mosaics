% Function within retinalMetamers
% By J. Freedland, 2020
%
% After calculating a neuron's receptive field (RF), a difference of gaussian
% (DoG) filter is created. This is used to normalize values in the
% naturalistic image's projection.
%
% INPUTS:   obj: structure from retinalMetamers
%               must contain:   obj.rfSigmaCenter (in microns)
%                               obj.rfSigmaSurround (in microns)
%                               obj.micronsPerPixel
%                               obj.videoSize (desired video size, in ARCMIN) 
%
% OUTPUTS:  RFFilter: receptive field filter
%%%

function [RFFilter] = calculateFilter(obj)

    % Convert neuron's RF to DOVES VH units.
    centerSigma_arcmin = cellMosaics.utils.changeUnits(obj.centerSigma,obj.micronsPerPixel,'um2arcmin');
    surroundSigma_arcmin = cellMosaics.utils.changeUnits(obj.surroundSigma,obj.micronsPerPixel,'um2arcmin');

    % Generate 2D gaussians
    centerGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],centerSigma_arcmin);
    surroundGaus = fspecial('gaussian',[obj.videoSize(1) obj.videoSize(2)],surroundSigma_arcmin);

    % Calculate difference of gaussians
    diffGaussian = centerGaus - surroundGaus;
    RFFilter = diffGaussian ./ max(diffGaussian(:)); % Normalize filter
end