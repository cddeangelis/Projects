clear;
%% Multi-layer Monte Carlo Simulation of NIR Light Propagation in the brain
% continuous wave functional near-infrared spectroscopy requires knowledge
% of photon paths in order to make measurements. It has no knowledge about
% the paths, thus Monte Carlo simulations are employed to make estimations.
% 
% this program simulates near-infrared light propagation through the brain. It is a
% simplified model with 5 planar layers, each with different properties
% corresponding to the scalp, skull, cerebral spinal fluid, gray matter,
% and white matter. This is intended to emulate functional near-infrared
% spectroscopy, where detectors and emitters are placed on the scalp. The
% program counts the photons that reach the detector and reports the
% mean path length of those photons.

% reference: Wang et al. "Monte Carlo modeling of light transport in
% multi-layered tissues" (1995)

% Limitations:
% no ray tracing
% single detector and emitter
% planar model
% identical refractive indicies (which is close to truth)

% detector is at (25,25) radius = 20
% SCALE is mm

% number of photons to be used
N = 1e7;

fprintf("For %0.0e photons....\n", N)

% layer depth
Tsc = 7.25; % scalp
Tsk = Tsc + 4; % skull
Tcb = Tsk + 2.73; % cerebral spinal fluid
Tgm = Tcb + 3.29; % gray matter
% white matter below Tgm

% a russian roulette is implemented to mimic absorption, more later
Wth = 0.0001; % weight threshold
m = 10; % russian roulette odds

% [scalp skull cbf gray white]
muA = [0.019, 0.019, 0.004, 0.02, 0.08]; % absorption coefficients
muS = [7.8, 7.8, 0.009, 9.0, 40.9]; % scattering coeffients 
g = [0.89, 0.89, 0.89, 0.89, 0.84]; % ansiotropy 
re = 1.37; % refractive index


paths = 0;
absorbed = 0;
escaped = 0;


for i = 1:N

    % initialization

    % cartesian coordiantes of photon
    x = 0;
    y = 0;
    z = 0;

    % direction cosines of photon
    ux = 0;
    uy = 0;
    uz = 1;

    % inital weight
    W = 1;
    
    pathlen = 0;
    currentlayer = 1;

    while 1

       % new step
        muI = muA(currentlayer) + muS(currentlayer); % interaction coefficient
        s = -log(rand)/muI; % probability distribution function  for step size
        pathlen = pathlen + s;
        
        % update coordiantes
        x = x + ux*s;
        y = y + uy*s;
        z = z + uz*s;

        dw = (muA(currentlayer)/muI) / W; % change in weight
        W = W - dw; % new weight

        % random deflection angle according to Henyey-Greenstein equation
        costheta = (1/(2*g(currentlayer)))*(1 + g(currentlayer)^2 - ((1-g(currentlayer)^2)/(1-g(currentlayer)+2*g(currentlayer)*rand))^2);
        sintheta = sqrt(1 - costheta^2);

        % random azimuthal angle, uniformly distributed
        fi = 2*pi*rand;


        % update angles, special equation if it is close to vertical
        if uz > 0.9999
            ux = sintheta*cos(fi);
            uy = sintheta*sin(fi);
            uz = costheta;
        elseif uz < -0.9999
            ux = sintheta*cos(fi);
            uy = -sintheta*sin(fi);
            uz = -costheta;
        else
            ux = sintheta*(ux*uz*cos(fi) - uy*sin(fi))/sqrt(1-uz^2) + ux*costheta;
            uy = sintheta*(uy*uz*cos(fi) + ux*sin(fi))/sqrt(1-uz^2) + uy*costheta;
            uz = -sqrt(1-uz^2)*sintheta*cos(fi) + uz*costheta;
        end

        % did it hit the detector?
        if z < 0
            if (x - 25)^2 + (y-25)^2 < 400
                tmp = size(paths);
                paths(tmp(2) + 1) = pathlen;
                break
            end
        end

        % did it escape?
        if z < 0
            escaped = escaped + 1;
            break
        end

        % update layer, which in turn updates optical properties
        if z > Tgm
            currentlayer = 5;
        elseif z > Tcb
            currentlayer = 4;
        elseif z > Tsk
            currentlayer = 3;
        elseif z > Tsc
            currentlayer = 2;
        else
            currentlayer = 1;
        end

        % Absorption
        % photon termination must follow conservation of energy without
        % skewing distribution of photon deposition, thus a Russian
        % roulette is used. 
        if W < Wth
            if rand < 1/m
                W = m*W;
            else
                absorbed = absorbed + 1;
                break % kill photon
            end
        end

        % if the photon didn't escape or get absorbed, it loops and takes
        % another step
    end

end

paths = paths(2:end);

tmp = size(paths);
nphotons = tmp(2);
mpl = mean(paths);

fprintf("\n%i photons reached the detector\n", nphotons)
fprintf("The mean path length was %f\n", mpl)
fprintf("%i photons were absorbed\n", absorbed)
fprintf("%i photons were lost\n", escaped);

