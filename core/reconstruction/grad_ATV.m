%Computes anisotropic total variation (ATV) and gradient of the ATV for an image,
%f. Assumes f is specified as a 2D image and returns the gradient as a 2D
%image as well. My implementation of ATV based on "A limited-angle CT reconstruction method based on
%anisotropic TV minimization" by Chen, Jin, Li and Wang. https://doi.org/10.1088/0031-9155/58/7/2119

%thetas are the angles for which data were acquired (used to weight directions), while alphas are
%directions used to calculate the ATV. See above paper for more detail.
function [TV,dTV] = grad_ATV(f, thetas, alphas)
    epsilon = 1e-6;

    [m, n] = size(f);
    TV = 0;
    dTV = zeros(size(f));
    omegas = zeros(1, numel(alphas));
    
    %computation of weighting parameters (see Section 2.4 of paper)
    for alpha_i = 1:numel(alphas)
        alpha = alphas(alpha_i);
        omegas(alpha_i) = compute_omega(thetas, alpha);
    end
    
    % See equations 17, 18
    omegas = omegas/norm(omegas, 2);
    
    %for each direction
    for alpha_i = 1:numel(alphas)
        alpha = alphas(alpha_i);
        omega = omegas(alpha_i);
        
        % Given alpha, calculate the unit vector e_alpha
        e_alpha = [cos(alpha/180*pi), sin(alpha/180*pi)]; 
        e_alpha_x = e_alpha(1); e_alpha_y = e_alpha(2); 
        
        ind_m0 = [m 1:m-1]; ind_n0 = [n 1:n-1];
        ind_m2 = [2:m 1]; ind_n2 = [2:n 1]; %assumes periodic BCs which should be ok if support is compact        
        ind_m1 = 1:m; ind_n1 = 1:n;
        
        % m,n
        diffx_1 = f(ind_m2, ind_n1) - f(ind_m1, ind_n1);
        diffy_1 = f(ind_m1, ind_n2) - f(ind_m1, ind_n1);
        
        % m-1, n
        diffx_0 = f(ind_m1, ind_n1) - f(ind_m0, ind_n1);
        diffy_0 = f(ind_m0, ind_n2) - f(ind_m0, ind_n1);        
        
        % m, n-1
        diffx_2 = f(ind_m2, ind_n0) - f(ind_m1, ind_n0);
        diffy_2 = f(ind_m1, ind_n1) - f(ind_m1, ind_n0);
        
        %dot product
        diff_1 = e_alpha_x*diffx_1 + e_alpha_y*diffy_1 + epsilon;
        diff_0 = e_alpha_x*diffx_0 + e_alpha_y*diffy_0 + epsilon;
        diff_2 = e_alpha_x*diffx_2 + e_alpha_y*diffy_2 + epsilon;        
        
        %norm of the above
        diffttl_1 = abs(diff_1);
        diffttl_0 = abs(diff_0);
        diffttl_2 = abs(diff_2);
        
        
        TV = TV + (omega * sum(sum(diffttl_1)));
        
        %gradient calculation
        alpha_i_dTV = ( (1./diffttl_1) .* (diff_1) * (-e_alpha_x - e_alpha_y) ) + ...
                      ( (1./diffttl_0) .* (diff_0) * e_alpha_x ) + ...
                      ( (1./diffttl_2) .* (diff_2) * e_alpha_y );


        dTV = dTV + (omega * alpha_i_dTV);        
    end    
end
