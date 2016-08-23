function rate = FF(adapt, target, settings)
fs = settings.sample_frequency;
eps = settings.eps; %0.125 s
minlen = settings.minlen; %0.5 s

tp = TP(adapt, target, eps*fs);
pv = FP(adapt, target, minlen*fs);
fn = length(target) - tp;
rate = 0;
if tp>0
    rate = 2*tp/(2*tp + pv + fn);
end


% Computation of TP (true positive)
function tp = TP(adapt, target, eps)
tp = 0;
for t = target
    if length(adapt(abs(adapt - t) <= eps))>0
        tp = tp + 1;      
    end
end

function fp = FP(adapt, target, minlen)
fp = 0;
for t = target
    idx = find(adapt(abs(adapt - t) <= minlen));
    adapt(idx) = [];  
end
fp = sum( abs(adapt(2:end) - adapt(1:end-1)) < minlen);

