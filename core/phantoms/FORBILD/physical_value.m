% function [ physphantom ] = physical_value(phantom, muW, muB)
% Converts relative attenuation values in the analytic phantom to physical
% mu-values for a desired energy.
% Input: phantom: the output from analytical_phantom.
%        muW: The attenuation coefficient for water at the desired energy
%        muB: The attenuation coefficient for bone at the desired energy
% Output:   physphantom is identical to the input phantom, except that the
%           sixth column of E (which consisted of attenuation coefficients relative to
%           water in the analytical phantom) now consists of interpolated, "true" attenuation
%           coefficients.

function [ physphantom ] = physical_value(phantom,muW,muB)
    physphantom.E=phantom.E;
    physphantom.C=phantom.C;
    nrows=size(phantom.E,1);
    shift=0;
    if (nrows > 71) shift=80;end
    if (nrows >= 97) physphantom.E(1:80,6)=muB-1.05*muW;end
    if (nrows == 71 || nrows == 151)
        physphantom.E(18+shift,6)=muB-1.05*muW;
        physphantom.E((19:71)+shift,6)=-muB;
    end

    physphantom.E(5+shift,6)=muB;
    physphantom.E(17+shift,6)=1.05*muW-muB;
    physphantom.E(6+shift,6)=-1.05*muW;
    physphantom.E(14+shift,6)=muB;
    physphantom.E([7 8 9 10 13 15 16]+shift,6)=muB-1.05*muW;
    j=[1 2 3 4 11 12];
    physphantom.E(j+shift,6)=muW*phantom.E(j+shift,6);
end