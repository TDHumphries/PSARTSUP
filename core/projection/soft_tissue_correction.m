%Apply soft tissue correction to measured data
%Input:
%sino: sinogram (pre-log measurements)
%energies: vector of energies contained in spectrum
%spectrum: beam spectrum values (column vector)
%mu_soft: soft tissue mu values at all energies. should be same size as spectrum.
%Eref: reference energy
%Output: simulated monoenergetic data (pre-log)
function m = soft_tissue_correction(sino, energies,spectrum, mu_soft, Eref )
 
    p = -log(sino);
    [m,n] = size(p);  
    l = zeros(m,n);
    mu_soft = mu_soft(:); spectrum = spectrum(:); %ensure column vectors
    ind = find(energies == Eref);
    mu_soft_ref = mu_soft(ind);    
    fprintf('Soft tissue correcting...\n');
    for i = 1:m
         fprintf('%d ',i); if( ~mod(i,20)) fprintf('\n'); end
        for j = 1:n
            func = @(t) p(i,j) + log((spectrum')*exp(-mu_soft*t));
            l(i,j) = fzero(func,p(i,j)/mu_soft_ref);
        end
    end
    fprintf('\n');
    m = exp(-l*mu_soft_ref);    %simulated monoenergetic data
    
end