%%% Code 2017, 2019 by Peter kasson
%%% Published in "Influenza hemifusion phenotype depends on membrane context: 
%%% differences in cell-cell and virus-cell fusion" by Zawada, Okamoto, and Kasson.
%%% Journal of Molecular Biology, 2018.
%%%
%%% Modified in "Detecting and controlling dye effects in single-virus fusion 
%%% experiments" by Rawle, Villamil Giraldo, Boxer, and Kasson.  Biophysical
%%% Journal, 2019.
%%%
%%% Extension and refactoring of code by Tijana Ivanovic:
%%% Ivanovic et al. "Influenza-virus membrane fusion: cooperative fold-back of
%%% stochastically induced hemaggutinin intermediates"
%%% https://doi.org/10.7554/eLife.00333
%%%
%%% Approach:
%%%  For each virion simulated:
%%%    Calculate a contact patch and compute an activation lag time for
%%%    each HA in the patch sampled randomly from an exponentially decay.
%%%    Also calculate an inactivation time at which point the HA inactivates.
%%%    Once a HA inactivates, it cannot contribute to hemifusion.
%%%
%%%    Sort particles by their lag times and calculate neighbors.
%%%    Calculate virion arrest time as activation of Narrest trimers.
%%%    Calculate hemifusion time as time when Nhemifusion neighboring 
%%%    trimers activate.
%%%    Calculate fusion time as time when Nfusion trimers activate and
%%%    undergo hemifusion->fusion transition, sampled from another
%%%    exponential decay.
%%%    Return vector of arrest, hemifusion, fusion times.

function [allarrest, allhemi, allfuse] = fusion_vs_inactivate(k_inact, ...
    k_hf, Nvirions, Ntrimers, Narrest, Nhemifusion, Finactive, Nfusion)

% now a parameter: k_hf = 0.0025;
Ndeadvirions = 0;
HA_spacing = 1;
kfuse = 0.55;  %from Floyd et al., 2008

% circular patch of Ntrimers elements arranged in a hexagonal lattice
[p,Nactual,a] = generate_patch(Ntrimers, HA_spacing);

% should each virion fuse, this is the waiting time after hemifusion
twait = random('exp', 1/kfuse, Nvirions, 1);

% Pre-allocate
allarrest(1:Nvirions) = NaN;
allhemi(1:Nvirions) = NaN;
allfuse(1:Nvirions) = NaN;


for ii = 1:Nvirions
    NotDeadFlag=0;
    tFusion = NaN;
    fl = 0;
    % sample randomly from an exponentially decaying funcion with rate k_hf;
    % get lag times for particles in order from 1st to Nactual-th
    p.flipped = random('exp', 1/k_hf, Nactual,1);
    p.lifetime = random('exp', 1/k_inact, Nactual,1);
    
    % sort data according to lag times
    [y,idx] = sortrows([p.id p.xy p.flipped],4); 
    p.id = y(:,1); p.xy = y(:,2:3); p.flipped = y(:,4);
    p.neighbors = p.neighbors(idx);   % Re-sort cell array.
    for kk = 1:length(p.neighbors)
        for jj = 1:length(p.neighbors{kk})
            newidx = find(p.neighbors{kk}(jj) == idx);
            p.neighbors{kk}(jj) = newidx;  
        end
    end
    
    p.inactive = zeros(Nactual,1);
    % inactivative according to lifetime calculations:
    for n=1:Nactual
        if p.flipped(n) > p.lifetime(n)
            % if the protein would flip after inactivation
            p.inactive(n) = 1;
        end
    end

    if Finactive==0 && k_inact <= 0.0
        tArrest = p.flipped(Narrest);
    else
        tArrest = NaN;
        r = randperm(Nactual);
        for jj = 1:ceil((Finactive * (Nactual)))
            p.inactive(r(jj))=1; 
        end
        
        for jj = 1: Nactual
            if p.inactive(jj) == 0
                fl = fl+1;
            end
            
            if fl==Narrest
                tArrest = p.flipped(jj);
                break
            end
        end       
    end
        
    for n = Nhemifusion:Nactual
        % after n flipped are there Nhemifusion neighbors?
        if p.inactive(n)==0
            if isaN2tuplet(p, n, Nhemifusion, a)
                NotDeadFlag = 1;  % or assign to random inactivation time?
                tHemifusion = p.flipped(n);
                break;
            end
        end   
    end

    for n = max(Nfusion, Nhemifusion):Nactual
        % after n flipped are there Nfusion neighbors?
        if p.inactive(n)==0
            if isaN2tuplet(p, n, Nfusion, a)
                % test p.flipped(n) > tinactive(n) or > NotDead if assign
                % above?
                NotDeadFlag = 1;
                tFusion = p.flipped(n) + twait(ii);
                break;
            end
        end   
    end    
    
    allarrest(ii) = tArrest;
    if NotDeadFlag == 0
        Ndeadvirions = Ndeadvirions+1;
        allhemi(ii) = NaN;
    else
        dt = tHemifusion-tArrest;
        allhemi(ii) = dt;
    end
    allfuse(ii) = tFusion;
end



    
