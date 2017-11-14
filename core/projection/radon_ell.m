%------- to generate the line integral in the ellipses--------
 function [RF] = radon_ell(E,theta,s)
%  This function computes the Radon transform of ellipses
% Input: E  as for MATLAB function phantom
%        theta as for MATLAB function radon
%        s: vector with radial coordinates corresponding to each row of RF
%        RF = line integral of the ellipses.
%

rho = E(:,1); u = E(:,2); v = E(:,3);
x = E(:,4); y = E(:,5); alpha = E(:,6)*pi/180;  %convert to radians
ne = length(rho);
RF = zeros(length(s),length(theta));
for j = 1:length(theta);
    phi = pi*theta(j)/180;
    omega = [cos(phi);sin(phi)];
    tmp =zeros(1,length(s));
    for mu = 1:ne
        a = (u(mu)*cos(phi-alpha(mu)))^2+(v(mu)*sin(phi-alpha(mu)))^2;
        test = a-(s-[x(mu);y(mu)]'*omega).^2;
        ind = test>0;
        tmp(ind) = tmp(ind)+rho(mu)*(2*u(mu)*v(mu)*sqrt(test(ind)))/a;
    end
    RF(:,j) = tmp.';
end
 end