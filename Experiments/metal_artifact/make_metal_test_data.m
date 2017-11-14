%% make_metal_test_data.m
% Code to generate data used in "Superiorized polyenergetic reconstruction algorithm 
% for reduction of metal artifacts in CT images" by Humphries & Gibali (http://faculty.uwashington.edu/thumphri/HG17.pdf)
%
% Contact T. Humphries (thumphri@uw.edu) with questions


addpath(genpath('../../core'));
addpath('../../../Libraries/fessler/irt');

setup;

%% 
%Phantom parameters
Ellpar =      [    0.9000    0.9000         0         0         0; %background
                    0.1000    0.1500   -0.4500   -0.4500         0; %bone 1
                    0.0750    0.0750   -0.4500    0.4500         0; %metal 1
                    0.1500    0.1000    0.4500    0.4500         0; %bone 1
                    0.0750    0.0750    0.4500   -0.4500         0; %metal 2
                    0.0460    0.0230   -0.0800    0.06         0;    %low contrast central features
                    0.0230    0.0230         0    0.06         0;
                    0.0230    0.0460    0.0600    0.06         0 ;            
                    0.0460    0.0230    0.0800    -0.06         0;
                    0.0230    0.0230         0    -0.06         0;
                    0.0230    0.0460    -0.0600    -0.06         0];                

%Spectrum and att data
kVp = 130;      
xrs=xray_read_spectra(sprintf('poly1,[%d]',kVp)); % use for I values
spectrum = xrs.sp/sum(xrs.sp); energies = xrs.en;   %normalize spectrum
mass_atten=xray_read_atten({'soft','bone','titanium'},1:kVp);
mass_dens = xray_read_dens({'soft','bone','titanium'});
soft_vals = mass_atten(:,1)'*mass_dens(1); %multiply by density; 
bone_vals = mass_atten(:,2)'*mass_dens(2);
titanium_vals = mass_atten(:,3)'*mass_dens(3); 


density = [soft_vals;   % intensity of elipses
    bone_vals-soft_vals; 
    titanium_vals-soft_vals; % attentuation coefficients
    bone_vals-soft_vals; 
    titanium_vals-soft_vals
    0.05*soft_vals
    0.05*soft_vals
    0.05*soft_vals
    -0.05*soft_vals
    -0.05*soft_vals
    -0.05*soft_vals];
%density = soft_vals;

theta = (-360:359)/4; % generate data on -90 to 90 degree arc
numpix = 400; dx = 30/numpix; % dx = pixelsize
Eref = 70;
s2=[-284.375:0.25:284.375]/numpix*2;
win = [0.18 0.23];

%% generate sinograms
 E = [density(:,Eref) Ellpar]; %Phantom at reference energy of 70kEV
 tmp = radon_ell(E,theta,s2);
 tmp = (tmp(1:4:end,:)+tmp(2:4:end,:)+tmp(3:4:end,:)+tmp(4:4:end,:))*0.25;
 sino_mono = exp(-0.5*dx*numpix*tmp);
 
 img_ref = phantom(E,numpix);
 
 sino_poly = zeros(size(sino_mono));
 fprintf('Generating poly data...');
 for j = 1:numel(spectrum)
     fprintf('%d ',j);
     E = [density(:,energies(j)) Ellpar];
     tmp = radon_ell(E,theta,s2);
     tmp = (tmp(1:4:end,:)+tmp(2:4:end,:)+tmp(3:4:end,:)+tmp(4:4:end,:))*0.25;
     sino_poly = sino_poly + spectrum(j)*exp(-0.5*dx*numpix*tmp);
 end
fprintf('\n');
 
%put data into correct format for reconstruction code. In particular,
%positive theta direction is CLOCKWISE and 0 degrees corresponds to taking
%measurements along x-axis (not y-axis) since theta refers to the view and
%not the angle that parameterizes the radon transform.
sino_mono = fliplr(sino_mono); sino_poly = fliplr(sino_poly);
theta = theta-90;

%soft tissue corrected data
 sino_poly_wc = soft_tissue_correction(sino_poly,energies,spectrum,soft_vals,Eref);
 
 save(sprintf('metal_test_data_%d.mat',kVp), 'sino_*','counts','Eref','numpix','dx','theta','img_ref', 'soft_vals');
 