% Function that creates monoenergetic sinogram data of the FORBILD phantom using the
% analytic projection operation. Input: FORBILD phantom object generated
% from analytical_phantom.m, attenuation value for water at that energy,
% attenuation value for bone at that energy, affine parameter s, angular
% parameter theta.
function sino =  project_mono_FORBILD(phantom, mu_w, mu_b, s, theta)
    P_phys = physical_value(phantom,mu_w,mu_b);
    
    num_s = numel(s);
    num_theta = numel(theta);
    scoord = ones(num_theta,1)*s;
    theta = transpose(theta)*ones(1,num_s);

    sino = exp(-line_integrals(scoord,theta,P_phys));
    sino = flipud(sino)';        %get it into the same format as for radon


end
