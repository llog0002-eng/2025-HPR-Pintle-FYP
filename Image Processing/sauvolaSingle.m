function BW = sauvolaSingle(I, window, k)
% SAUVOLASINGLE - Single-stage Sauvola binarization
% I : grayscale image (normalized or raw)
% window : [M N] neighborhood (default [3 3])
% k : Sauvola sensitivity (default 0.34)

if nargin < 3
    k = 0.34;
end
if nargin < 2
    window = [3 3];
end

if ~ismatrix(I)
    error('Input image must be a 2D array.');
end
I = double(I);

% Compute local mean and standard deviation using averagefilter
m  = averagefilter(I, window, 'replicate');
m2 = averagefilter(I.^2, window, 'replicate');
s  = (m2-m.^2).^0.5;

% Sauvola threshold
R = max(s(:));
T = m .* (1 + k * (s./R - 1));

% Apply threshold
BW = (I > T);
end