%computation of weighting parameter omega used in calculation of ATV. See
%Section 2.4 of "A limited-angle CT reconstruction method based on
%anisotropic TV minimization" by Chen, Jin, Li and Wang. https://doi.org/10.1088/0031-9155/58/7/2119
%Input:
%   thetas: angles corresponding to the projection data being used
%   alpha: direction in which aTV is being calculated
%Output
%   omega: the corresponding weight
function [omega] = compute_omega(thetas, alpha)
    omega = 0;
    
    for theta_j = 1:numel(thetas)
        theta = thetas(theta_j);
  
        theta = mod(theta,360); alpha = mod(alpha,360);
        if theta > 180
            theta = theta - 180;
        end
        if alpha > 180
            alpha = alpha - 180;
        end

        br_theta = theta;% + 90; changed this for MIRT angular convention
        if (br_theta > 180)
            br_theta = br_theta - 180;
        end
        
        omega = omega + abs(sin(alpha*pi/180 - br_theta*pi/180));
    end    
end