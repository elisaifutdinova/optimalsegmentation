 function[adapt_borders] = segmentation(data, settings)

window_length = settings.window_lenght;
step = settings.window_step;
fs = settings.sample_frequency;
coef = settings.threshold_coef;
trsh_window = settings.trsh_window;


% moving windows segemntation algorithm
window = ceil( window_length * fs); % (samples)
window_subnum = ceil(window/step); % (samples)

% buffering
data_buffered = buffer(data, window, [window-step], 'nodelay');
total_segmments = size(data_buffered, 2);

% functions
func_std = varri(data_buffered, settings);
clear data_buffered;

% indexes for two windows
first_window_indexes = 1 : total_segmments-window_subnum;
second_window_indexes = first_window_indexes + window_subnum;

% difference between two windows
func_std_diff = abs(func_std(first_window_indexes) - func_std(second_window_indexes));
func_std_diff = func_std_diff/max(func_std_diff);
MPH = std(func_std_diff);

h = hann(trsh_window*step);
h = h/sum(h);
threshold = coef*conv(func_std_diff, h, 'same');
idx = find(func_std_diff < threshold);
func_std_diff(idx) = threshold(idx);%threshold(idx);


MPD = window/100; % (s)
MPH = std(func_std_diff);
MPD = floor((MPD*fs) / step);
MPD = ceil(window/step);
[~, locs] = findpeaks(func_std_diff, 'MINPEAKDISTANCE', MPD, 'MINPEAKHEIGHT', MPH);
adapt_borders = window + (locs-1)*step;
 

function values = varri(data, options)
k1 = options.varr_parameter1;
k2 = options.varr_parameter2;
data = data';
for k = 1:length(data(:, 1))
    x = data(k, :);
    values(k) = k1*sum(abs(x)) + k2*sum(abs(x - [x(2:end) 0]));
end



    
    
    
    
    

