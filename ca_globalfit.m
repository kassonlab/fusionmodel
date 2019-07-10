%%% Code 2019 by Peter Kasson
%%% Calculate fit parameters for cellular-automaton fusion model
%%% against a series of experiments and hemifusion times.
%%% Published in Rawle, Villamil Giraldo, Boxer, Kasson,
%%% Biophysical Journal, 2019.
%%% Underlying model from Ivanovic et al., 2013; modified by
%%% Zawada, Okamoto, Kasson, 2018.


function [nll, params] = ca_globalfit(startvals, params_to_vary, ...
    reactants, dependent_step, hf_times, efflist)
    % fit parameters for cellular automaton model against a set of
    % hemifusion times.  Return negative log-likelihood and best params.
    gs = GlobalSearch;
    gs.NumTrialPoints = 5000;
    % can set options
    obj = @(x)nll_par_sub(x, startvals, params_to_vary, reactants, ...
        dependent_step, hf_times, efflist);
    problem = createOptimProblem('fmincon', 'x0', ...
        startvals(params_to_vary), 'objective', obj, 'lb', ...
        0*startvals(params_to_vary), 'ub', 100*startvals(params_to_vary));
    [params, nll] = run(gs, problem);
    % [params, nll, ~, ~] = particleswarm(obj, 2);
end

function nll=nll_par_sub(x, static_param, var_idx, reactants, ...
    dependent_step, hf_times, efflist)
    static_param(var_idx) = x;
    nll = global_ca_nll(static_param, reactants, dependent_step, ...
        hf_times, efflist);
end

function nll=nll_ca(params, reactants, dependent_step, hf_times, eff)
    % calculate cellular-automaton-model NLL
    params(dependent_step) = params(dependent_step) * reactants;
    cellparam = num2cell(params);
    [~, hemi, ~] = fusion_vs_inactivate(cellparam{:});
    if isnan(hemi)
        % all unfused; do efficiency only
        nll = -1 * log(eff) * length(hemi);
    else
        [~, likelihoods] = ksdensity(hemi, hf_times);
        nll = -1 * sum(log(likelihoods(isfinite(likelihoods))));
        nll = nll - log(eff) * nnz(isnan(likelihoods));
    end
end

function nll=global_ca_nll(params, reactant_list, dependent_step, hf_timelist, efflist)
    % calculate NLL across multiple cellular automata
    nll = 0;
    for i=1:size(hf_timelist, 1)
        nll = nll + nll_ca(params, reactant_list(i), dependent_step, ...
            hf_timelist{i}, efflist(i));
    end
end

    
