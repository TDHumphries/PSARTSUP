%Computes total variation and gradient of the total variation for an image,
%f. Assumes f is specified as a 2D image and returns the gradient as a 2D
%image as well.
function [TV,dTV] = grad_TV(f)
epsilon = 1e-6;
[m,n] = size(f);

%create indices in x and y directions, offset by 1 in either direction
ind_m1 = 1:m; ind_n1 = 1:n;
ind_m2 = [2:m 1]; ind_n2 = [2:n 1];         %assumes periodic BCs which should be ok if support is compact
ind_m0 = [m 1:m-1]; ind_n0 = [n 1:n-1];
diff1 = (f(ind_m2,ind_n1)-f(ind_m1,ind_n1)).^2; %difference in y-direction
diff2 = (f(ind_m1,ind_n2)-f(ind_m1,ind_n1)).^2; %difference in x-direction
diffttl = sqrt(diff1+diff2+epsilon^2);          %gradient norm in every pixel
TV = sum(sum(abs(diffttl)));

dTV = -1./diffttl .* (f(ind_m2,ind_n1)-2*f(ind_m1,ind_n1)+f(ind_m1,ind_n2)) + ...
        1./diffttl(ind_m0,ind_n1) .* (f(ind_m1,ind_n1)-f(ind_m0,ind_n1)) + ...
        1./diffttl(ind_m1,ind_n0) .* (f(ind_m1,ind_n1)-f(ind_m1,ind_n0));       %gradient of TV
%dTV = dTV(:);
end