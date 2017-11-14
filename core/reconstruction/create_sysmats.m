function [A, Dinv, Minv] = create_sysmats(geom, arc, numbins, numangles, start_angle, numpix,dx,fb_params, ns, use_weights, weights)

    A = cell(ns,1); Dinv=cell(ns,1); Minv=cell(ns,1);
    ig = image_geom('nx',numpix,'dx',dx);
    da = arc/numangles; %angular delta
    
%     fname = sprintf('sysmat_%dpx_%dviews_%dsubsets.mat',numpix,numangles,ns);
%     if (exist(fname,'file'))
%         fprintf('Loading sysmat from file\n');
%         load(fname);
%         return
%     end
        
    
    parfor j = 1:ns
        fprintf('Creating sysmat %d of %d\n',j,ns);
        sa_sub = start_angle + da*(j-1);    %start angle for this subset of views    
        switch(geom)
            case 'fan'
                sg = sino_geom('fan','orbit_start',sa_sub,'orbit', arc, 'nb',numbins,'na',numangles/ns,'ds',fb_params.ds,'dso',fb_params.dso,'dod',fb_params.dod,'units','cm','strip_width','ds');
            case 'par'
                sg = sino_geom('par','orbit_start',sa_sub,'orbit', arc, 'nb',numbins,'na',numangles/ns,'dr',dx,'strip_width','d');
            otherwise
                error('Invalid geometry specified');
        end

        mat = Gtomo2_strip(sg,ig);
        
        %Compute equivalent to matrix D
        if use_weights 
            v = weights(:,j:ns:numangles); v = v(:) + eps;
        else
            v = ones(numbins*numangles/ns,1)+eps;
        end
        tmp = mat'*v;
        Dinv{j} = reshape(tmp,numpix,numpix);
        
        %compute equivalent to matrix M
        tmp = mat*ones(numpix*numpix,1)+eps; 
        Minv{j} = reshape(tmp,numbins,numangles/ns);
        A{j} = mat;
        %clear tmp mat       
    end
    %save(fname,'A','Dinv','Minv','-v7.3');

end