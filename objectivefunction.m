function f=objectivefunction(Swarm,data)

settings.sample_frequency = data.sample_frequency;
settings.eps = data.eps;
settings.minlen = data.minlen;

for i=1:size(Swarm,1),
    
    x=Swarm(i,:);
    
    %% settings
    settings.window_step = round(x(1)); % in samples // 5-100
    settings.window_lenght = x(2);      % in seconds // 0.05-10
    settings.threshold_coef = x(3);     % 0 - 2
    settings.varr_parameter1 = x(4);    % 0 - 100
    settings.varr_parameter2 = x(5);    % 0 - 100
    settings.trsh_window = x(6);        % 5-20
    
    l = length(data.data);
    cr = zeros(1, l);
    
    for j = 1:l
        settings.fs = data.data(j).fs;
        d = data.data(j);
        % adaptive segmentation
        adapt = segmentation(d.signal, settings);
        % criteria value
        cr(j) = FF(adapt, d.target, settings);
    end
   
    f(i) = -mean(cr);
    disp([x f(i)]);

end

