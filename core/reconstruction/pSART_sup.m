function [f,k,res,err] = pSART_sup(sino, geom_params, recon_params, poly_params, sup_params)

    lp = -log(sino);
    dims = size(sino);
    numbins = dims(1); numangles = dims(2);
    
    
    %Read in parameters
    assert(logical(exist('geom_params','var')),'Must specify geom_params!')
    [geom, arc, start_angle, numpix, dx, fb_params] = get_geom_params(geom_params);
    
    if ~exist('recon_params','var'), recon_params = []; end
    [max_its, ns, x_true, x0, beta, use_weights] = get_recon_params(recon_params, numangles,numpix);
        
    use_poly = false; use_iodine = false;
    if (exist('poly_params','var') && ~isempty(poly_params))
        use_poly=true;
        [energies,spectrum,Eref, material_att, material_names, mask, use_iodine] = get_poly_params(poly_params);
    end
    
    use_sup = false;
    if exist('sup_params','var') && ~isempty(sup_params)
        use_sup = true;
        [gamma, alpha_init, N, objfun, epsilon_target] = get_sup_params(sup_params);
        alpha = alpha_init; 
    end
    
    if use_weights
        weights = sqrt(sino);
    else
        weights = [];
    end
    
    %generate system matrices and normalization factors
    [A, Dinv, Minv] = create_sysmats(geom, arc, numbins, numangles, start_angle, numpix,dx,fb_params, ns, use_weights, weights); 
    
    res = []; err = []; 
    f = x0;
    print_params(use_poly,use_iodine,use_sup,use_weights);
    
    %begin iteration
    for k = 1:max_its; %max number of iterations
        fprintf('It #%d\n',k)
        
        %Superiorization loop
        if(use_sup)
            fprintf('\tTV loop:');
            [g,dg] = objfun(f);           
            g_old = g; dg = -dg/(norm(dg(:),2)+eps);
            for j = 1:N
                fprintf('%d',j);        
                while true
                    fprintf('.');  %used to indicate progress of TV minimization/how many time it has to backtrack
                    f_tmp = f + alpha*dg;
                    g_new = objfun(f_tmp);
                    alpha = alpha*gamma;
                    if g_new <= g_old
                        %TV_old = TV_new; %comment out for less rigorous condition
                        f = f_tmp;
                        break
                    end
                end    
                [~,dg] = objfun(f);
                dg = -dg/(norm(dg(:),2)+eps); %ensures no division by zero in constant vector case
            end
            fprintf('\n');
        end
        
        %SART step
       tic;
       for j = 1:ns
           index = j:ns:numangles; 
           dims2 = [numbins numangles/ns]; %size of subset sinogram
           if (use_poly)
                fp = project_poly(f,Eref,energies,spectrum,material_att,material_names,A{j},dims2,use_iodine,mask);
           else
                fp = exp(-A{j}*f(:)); fp = reshape(fp,dims2);
           end
           diffs = (lp(:,index) + log(fp));
           diffs = diffs./Minv{j};
           if use_weights
               diffs = weights(:, index).*diffs;
           end
           bp = (A{j}')*diffs(:);
           bp = reshape(bp,numpix,numpix);    
           f = f + beta*bp ./ Dinv{j};
           f = max(f,eps); %ensure positivity
        end
        err(k) = norm(f(:) - x_true(:),2);
        %compute total projection to compute epsilon value
        for j = 1:ns
           index = j:ns:numangles; 
           dims2 = [numbins numangles/ns]; %size of subset sinogram
           if (use_poly)
                fp(:,index) = project_poly(f,Eref,energies,spectrum,material_att,material_names,A{j},dims2,use_iodine,mask);
           else
                tmp = exp(-A{j}*f(:)); fp(:,index) = reshape(tmp,dims2);
           end
        end
        epsilon = norm(lp(:) + log(fp(:)),2);
        res(k) = epsilon;
            toc;
        fprintf('it %d: res %e \t err %e\n',k, res(k),err(k));
        if (use_sup) && (epsilon < epsilon_target)
            return
        end

    end
    
end