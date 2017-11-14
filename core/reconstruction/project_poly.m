function proj = project_poly(f,Eref,energies,spectrum,material_att,material_names,A,dims,use_iodine,mask)

    [M,~] = size(material_att);
    assert(~use_iodine || M >= 4,'Must have at least four materials when using iodine (e.g. air, soft tissue, iodine, bone)');
    i = (energies == Eref); base_materials_ref = material_att(:,i);
    vals = zeros([dims M]);
    
    %find index of iodine attenuation coefficient
    if(use_iodine)
        K = find(strncmp(material_names,'iodine',6));
        assert(~isempty(K),'Iodine not specified as one of the materials with use_iodine');
        assert(numel(K) == 1,'More than one element specified with iodine in name');
        assert(K >= 3,'Need two tissues below iodine (e.g. air and soft tissue');
    else
        K = -1; %ensures that iodine method will not be used
        mask = [];
    end
    
    for j = 1:M
        tmp = zeros(size(f));
        switch j
            case 1                                              %should be air
                %attenuation greater than this tissue
                ind = (f < base_materials_ref(2));      
                tmp(ind) = (base_materials_ref(2)-f(ind))/(base_materials_ref(2)-base_materials_ref(1));
            case M
                %attenuation greater than previous tissue in list. Depends
                %on whether previous material is iodine or not.
                if (M-1 == K)
                    ind = (f >= base_materials_ref(M-2)) & mask ;
                    tmp(ind) = (f(ind) - base_materials_ref(M-2))/(base_materials_ref(M)-base_materials_ref(M-2));
                else
                    ind = (f >= base_materials_ref(M-1));
                    tmp(ind) = (f(ind) - base_materials_ref(M-1))/(base_materials_ref(M)-base_materials_ref(M-1));
                end
            case K - 1                                 %tissue before iodine -- typically soft tissue
                %attenuation less than this tissue
                ind = (f < base_materials_ref(K-1));  
                tmp(ind) = (f(ind)-base_materials_ref(K-2))/(base_materials_ref(K-1)-base_materials_ref(K-2));
                
                %attenuation greater than this tissue -- iodine mix
                ind = (f < base_materials_ref(K) & f >= base_materials_ref(K-1))&(~mask);  %mask is 0 for soft/iodine regions
                tmp(ind) = (base_materials_ref(K)-f(ind))/(base_materials_ref(K)-base_materials_ref(K-1));  
                
                %attenuation greater than this tissue -- bone (?) mix
                ind = (f < base_materials_ref(K+1) & f >= base_materials_ref(K-1)) & mask;  %mask is 1 for soft/bone regions
                tmp(ind) = (base_materials_ref(K+1)-f(ind))/(base_materials_ref(K+1)-base_materials_ref(K-1));  
            case K
                %attenuation greater than previous tissue in list
                ind = (f >= base_materials_ref(K-1))&(~mask);  
                tmp(ind) = (f(ind) - base_materials_ref(K-1))/(base_materials_ref(K)-base_materials_ref(K-1));  
            otherwise
                 %attenuation greater than this tissue
                 ind = (f < base_materials_ref(j+1) & f >= base_materials_ref(j));       
                 tmp(ind) = (base_materials_ref(j+1)-f(ind))/(base_materials_ref(j+1)-base_materials_ref(j));
                 %attenuation less than this tissue
                 ind = (f < base_materials_ref(j) & f >= base_materials_ref(j-1));
                 tmp(ind) = (f(ind)-base_materials_ref(j-1))/(base_materials_ref(j)-base_materials_ref(j-1));
        end
          vals(:,:,j) = reshape(A*tmp(:),dims);
    end
    
    %Implementation of Equation 12
    proj = zeros(dims);
    for j = 1:numel(energies)
        base_materials_e = material_att(:,j);
        tmp = zeros(dims);
        for k = 1:M
            tmp = tmp + vals(:,:,k)*base_materials_e(k);
        end
        proj = proj + spectrum(j)*exp(-tmp); 
    end
end

