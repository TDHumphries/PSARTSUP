function [max_its, ns, x_true, x0, beta, use_weights] = get_recon_params(params,na,numpix)
    
    fprintf('Reading reconstruction parameters...\n');
    
    %maximum number of iterations
    if ~isfield(params,'max_its')
        fprintf('\tMax iterations not specified in recon_params.max_its -- using 20\n');
        max_its = 20;
    else
        max_its = params.max_its;
    end
    
    %number of subsets
    if ~isfield(params,'ns')
        fprintf('\tNumber of subsets not specified in recon_params.ns -- using 1\n');
       	ns = 1;
    else
        ns = params.ns;
        assert(mod(na,ns) == 0,'Number of subsets (%d) does not divide number of angles (%d)',ns,na);
    end    

    %true image
    if ~isfield(params,'x_true')
        fprintf('\ttrue image not specified in recon_params.x_true -- error will not be calculated correctly\n');
       	x_true = zeros(numpix);
    else
        x_true = params.x_true;
    end    
    
    %starting image
    if ~isfield(params,'x0')
        fprintf('\tInitial guess not specified in recon_params.x0 -- using zeros\n');
       	x_true = zeros(numpix);
    else
        x0 = params.x0;
    end
    
   if ~isfield(params,'beta')
        fprintf('\tRelaxation parameter not specified in recon_params.beta -- using 1\n');
       	beta = 1.0;
   else
        beta = params.beta;
        assert(beta < 2 && beta > 0, 'relaxation parameter beta = %3.2f not between 0 and 2',beta);
   end    

   if ~isfield(params,'use_weights')
        fprintf('\tWeighting not specified in recon_params.use_weights -- assuming no weighting\n');
       	use_weights = false;
    else
        use_weights = params.use_weights;
    end        
end