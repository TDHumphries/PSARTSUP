function [gamma, alpha_init, N, objfun, epsilon_target] = get_sup_params(params)
    
    fprintf('Reading superiorization parameters...\n');
    
    if ~isfield(params,'gamma')
        fprintf('\tStep size parameter not specified in sup_params.gamma -- using 0.999\n');
       	gamma = 0.999;
    else
        gamma = params.gamma;
    end    

    if ~isfield(params,'alpha_init')
        fprintf('\tInitial alpha parameter not specified in sup_params.gamma -- using 1\n');
       	alpha_init = 1;
    else
        alpha_init = params.alpha_init;
    end    
    
    if ~isfield(params,'N')
        fprintf('\tNumber of inner loop iterations not specified in sup_params.N -- using 20\n');
       	N = 20;
    else
        N = params.N;
    end    
    
  	assert(isfield(params,'objfun'),'Must specify objective function for superiorization');
    objfun = params.objfun;
    
   if ~isfield(params,'epsilon_target')
        fprintf('\tTarget epsilon not specified in sup_params.N -- using 0\n\tMethod will iterate until max_its is attained.');
       	epsilon_target = 0;
    else
        epsilon_target = params.epsilon_target;
    end    
    
    
end