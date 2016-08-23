function results=mm_pso(options, data)
GM = -Inf; % Global minimum (used in the stopping criterion)
ErrGoal = 1e-10; % Desired accuracy
% 

% Initializing variables
success = 0; % Success Flag
iter = 0;   % Iterations' counter
fevals = 0; % Function evaluations' counter

% Using params---
% Determine the value of weight change
w_start = options.w_start;   %Initial inertia weight's value
w_end = options.w_end;       %Final inertia weight
w_varyfor = floor(options.w_varyfor*options.Iterations); %Weight change step. Defines total number of iterations for which weight is changed.
w_now = w_start;
inertdec = (w_start-w_end)/w_varyfor; %Inertia weight's change per iteration

% Initialize Swarm and Velocity
SwarmSize = options.SwarmSize;
Swarm = rand(SwarmSize, options.Dim).*repmat(options.ub-options.lb,SwarmSize,1) + repmat(options.lb,SwarmSize,1);
VStep = rand(SwarmSize, options.Dim);
f2eval = options.f2eval; %The objective function to optimize.

%Find initial function values.
%keyboard;
fSwarm= feval(f2eval, Swarm, data);

% Initializing the Best positions matrix and
% the corresponding function values
PBest = Swarm;
fPBest = fSwarm;

% Finding best particle in initial population
[fGBest, g] = min(fSwarm);
lastbpf = fGBest;
Best = Swarm(g,:); %Used to keep track of the Best particle ever
fBest = fGBest;
history = [];

if options.Neighbor
    % Define social neighborhoods for all the particles
    for i = 1:SwarmSize
        lo = mod(i-options.Nhood:i+options.Nhood, SwarmSize);
        nhood(i,:) = [lo];
    end
    nhood(find(nhood==0)) = SwarmSize; %Replace zeros with the index of last particle.
end


%
%                  THE  PSO  LOOP                          %
%
while( (success == 0) & (iter <= options.Iterations) )
    iter = iter+1;
    
    % Update the value of the inertia weight w
    if (iter<=w_varyfor) & (iter > 1)
        w_now = w_now - inertdec; %Change inertia weight
    end
    
    % The PLAIN PSO %
    
    % Set GBest
    A = repmat(Swarm(g,:), SwarmSize, 1); %A = GBest. repmat(X, m, n) repeats the matrix X in m rows by n columns.
    B = A; %B wil be nBest (best neighbor) matrix
    
    % use neighborhood model
    % circular neighborhood is used
    if options.Neighbor
        for i = 1:SwarmSize
            [fNBest(i), nb(i)] = min(fSwarm( find(nhood(i)) ));
            B(i, :) = Swarm(nb(i), :);
        end
    end
        
    % Generate Random Numbers
    R1 = rand(SwarmSize, options.Dim);
    R2 = rand(SwarmSize, options.Dim);
    
    % Calculate Velocity
    if ~options.Neighbor %Normal
        VStep = w_now*VStep + options.c1*R1.*(PBest-Swarm) + options.c2*R2.*(A-Swarm);
    else %With neighborhood
        R3 = rand(SwarmSize, options.Dim); %random nos for neighborhood
        VStep = w_now*VStep + options.c1*R1.*(PBest-Swarm) + options.c2*R2.*(A-Swarm) + options.c3*R3.*(B-Swarm);
    end
    
    %Vmax must be a matrix, because Vmax is different for different
    %components:
    Vmax=repmat(options.Vmax,SwarmSize,1);
    
    % Apply Vmax Operator for v > Vmax
    changeRows = VStep > Vmax;
    VStep(changeRows) = Vmax(changeRows);
    % Apply Vmax Operator for v < -Vmax
    changeRows = VStep < -Vmax;
    VStep(changeRows) = -Vmax(changeRows);
    
    % ::UPDATE POSITIONS OF PARTICLES::
    Swarm = Swarm + options.Chi * VStep;    
    xlimitsl=repmat(options.xlimits(1,:),SwarmSize,1);
    xlimitsu=repmat(options.xlimits(2,:),SwarmSize,1);
    VStep(Swarm<xlimitsl)=-VStep(Swarm<xlimitsl);
    VStep(Swarm>xlimitsu)=-VStep(Swarm>xlimitsu);
    Swarm(Swarm<xlimitsl)=xlimitsl(Swarm<xlimitsl);
    Swarm(Swarm>xlimitsu)=xlimitsu(Swarm>xlimitsu);
    %Constraint handling for x(1)<512*x(2);
    fs=data.sample_frequency;
    S1=Swarm(:,1);
    S2=Swarm(:,2);
    iv=S1>fs*S2;
    Swarm(iv,1)=(fs^2*S1(iv)+fs*S2(iv))/(fs^2+1);
    Swarm(iv,2)=Swarm(iv,1)/fs;
    VStep(iv,[1 2])=-VStep(iv,[1 2]);
    
    
    
    subplot(1,2,1);
    cla;
    hold on;
    line([0 10],  [0 5120]);
    plot(Swarm(:,2),Swarm(:,1),'g+');    
    plot(Swarm(iv,2),Swarm(iv,1),'r+');    
    axis([0.05 10 5 100]);
    xlabel('x_2 - window_length');
    ylabel('x_1 - window step');
    drawnow;
    
    
    % Evaluate new Swarm
   fSwarm= feval(f2eval, Swarm, data);%+Penalty;

    
    % Updating the best position for each particle
    changeRows = fSwarm < fPBest;
    fPBest(changeRows) = fSwarm(changeRows);
    PBest(changeRows, :) = Swarm(changeRows, :);
    
   
    % Updating index g
    [fGBest, g] = min(fPBest);

    %Update Best. Only if fitness has improved.
    if fGBest < lastbpf
        [fBest, b] = min(fPBest);
        Best = PBest(b,:);
    
    
        

    end
    
    history.Swarm(:,:,iter)=Swarm;
    history.fSwarm(:,:,iter)=fSwarm;
    history.Best(:,:,iter)=Best;
    history.fBest(:,:,iter)=fBest;
    
    
    
    
    
    if (rem(iter, options.DispInt) == 0)
        disp('************************');
        disp(sprintf('%4d\t\t\t%.5g\t\t\t%5d', iter, fBest));
        
    end
    
   subplot(1,2,2);hold on;
   plot(iter, fSwarm,'k+');
   plot(iter, fBest,'ro');
   xlabel('Iterations');
   ylabel('Cost');
   drawnow;
    
    lastbpf = fGBest; %To be used to find Best

    %
    %Just for the case of a sudden crash of the run, save partial results
    %save result_partial_pso;    
end
    

[fxmin, b] = min(fPBest);
xmin = PBest(b, :);

results.xmin=xmin;
results.fmin=fxmin;
results.history=history;

