%% make_FORBILD_data.m
% Code to generate data used in "Superiorized algorithm for reconstruction
% of CT images from sparse-view and limited-angle polyenergetic data" by
% Humphries, Winn and Faridani. https://arxiv.org/abs/1701.03396
%
% Contact T. Humphries (thumphri@uw.edu) with questions


%% Load libraries
addpath(genpath('../../core/'));
addpath('../../../Libraries/fessler/irt'); %put your own path to IRT here

setup;

%% Phantom/projection Geometry

 %800 x 800
 numpix = 800; dnp = (numpix-1)/2;
 dx = 30/numpix;                    %pixel width for total size of 30 cm
 s = (-567.375:0.25:567.375)*dx; %generate 4 lines for every measurement
                                 %this set of s values is chosen for compatibility with Matlab's radon command.
                                 %radon generates 1135 measurements for 800x800 image, so this generates 4x as
                                 %many, equally spaced between -567.5 and +567.5, and scaled by the pixel width.
 
theta_rad = (0:2879)*2*pi/2880-pi;  %generate 2880 angles between -pi and pi
theta_deg = theta_rad*180/pi;       %convert to degrees

%arrays used by FORBILD's projection code
x = ((0:numpix-1)-dnp)*dx;
y = ((0:numpix-1)-dnp)*dx;
xcoord = ones(numpix,1)*x;
ycoord = transpose(y)*ones(1,numpix);


%% Spectral data

xrs=xray_read_spectra('poly1,[130]'); % spectrum values
spectrum = xrs.sp/sum(xrs.sp); energies = xrs.en; %normalize spectrum to have integral of 1
mass_atten=xray_read_atten({'soft','bone'},[1:130]);
soft_vals = mass_atten(1:130,1)*1.06; %multiply by density to get linear attenuation coefficient; 
bone_vals = mass_atten(1:130,2)*1.92;

Eref = 70; %reference energy

mu_bone = bone_vals(Eref);
mu_soft = soft_vals(Eref)/1.05;   %normalize so that soft tissue val is the main background in FORBILD phantom
win = [0.195 0.215];

%% Generate monoenergetic data. Since we use analytic line integrals, this produces data slightly inconsistent with the MIRT system matrix.

 P = analytical_phantom(0,1); %FORBILD phantom with right ear insert
 P_phys =  physical_value(P, mu_soft ,mu_bone);
 img_70kev = discrete_phantom(xcoord,ycoord,P_phys); %generate discrete phantom at reference energy

%monoenergetic data. Average over four lines 
sino_mono_analytic = project_mono_FORBILD(P,mu_soft,mu_bone,s,theta_rad);
sino_mono_analytic = exp(- 0.25*(-log(sino_mono_analytic(1:4:end,:)) - log(sino_mono_analytic(2:4:end,:))-log(sino_mono_analytic(3:4:end,:))-log(sino_mono_analytic(4:4:end,:))));


recon_mono = iradon(-log(sino_mono_analytic(:,1:1440)),theta_deg(1:1440),'Hamming',numpix)/dx;
 
figure(1); 
subplot(1,2,1); imagesc(img_70kev,win); axis image; colormap gray;
subplot(1,2,2); imagesc(recon_mono,win); axis image; colormap gray;
 
%% Generate polyenergetic data. This takes a while. 
 
sino_poly_analytic = zeros(size(sino_mono_analytic));
  for j = 1:numel(energies) 
      disp(j);
      E = energies(j);
      mu_soft = soft_vals(E)/1.05; mu_bone = bone_vals(E);
      sino_tmp = project_mono_FORBILD(P,mu_soft,mu_bone,s,theta_rad);
      sino_tmp = exp(- 0.25*(-log(sino_tmp(1:4:end,:)) - log(sino_tmp(2:4:end,:))-log(sino_tmp(3:4:end,:))-log(sino_tmp(4:4:end,:))));
      sino_poly_analytic = sino_poly_analytic + spectrum(j)*sino_tmp;
  end
% 

%% Convert mono and poly data to convention used by Michigan Image Reconstruction toolbox
% Code was originally implemented for compatibility with Matlab's radon
% command. However, MIRT uses different angular convention. Specifically,
% positive rotation direction is clockwise (not counterclockwise) and zero
% degrees corresponds to taking line integrals in the positive x-direction,
% not the positive y-direction as in radon.

sino_mono_analytic = fliplr(sino_mono_analytic); sino_poly_analytic = fliplr(sino_poly_analytic); %reverses direction of rotation
sino_mono_analytic = sino_mono_analytic(:,[721:2880 1:720]); sino_poly_analytic = sino_poly_analytic(:,[721:2880 1:720]);

 save FORBILD_data img_70kev sino_mono_analytic sino_poly_analytic theta_deg numpix dx spectrum energies soft_vals bone_vals
 
%% Add poisson noise to data simulating the count rate indicated below

load FORBILD_data

counts = 4e6;
sino_mono_noisy = 1e12/counts*imnoise(counts/1e12*sino_mono_analytic,'poisson');
sino_poly_noisy = 1e12/counts*imnoise(counts/1e12*sino_poly_analytic,'poisson');

save FORBILD_data img_70kev sino_mono_analytic sino_poly_analytic theta_deg numpix dx spectrum energies soft_vals bone_vals sino_mono_noisy sino_poly_noisy

%% generate spectral approximation used for inconsistent data reconstructions
load FORBILD_data

figure(1); clf; plot(energies,spectrum,'b-'); hold on;
selectedenergies = [3 5 8 10 11 12 16 20 24 28 32 40 48 56:60 66:70 100 130]; %specifically chosen energies from the spectrum
energies_r = energies(selectedenergies); spectrum_r = spectrum(selectedenergies);
plot(energies_r,spectrum_r,'r-x'); hold off;
diffs = diff(energies_r);
spectrum_r = ([ spectrum_r(1:end-1).*diffs; 0] + [0; spectrum_r(2:end).*diffs])/2; %applies composite trapezoid rule
spectrum_r = spectrum_r/sum(spectrum_r); %renormalize to 1
spectrum_r = spectrum_r(2:end-1);       %end values of spectrum are zero and can be omitted
energies_r = energies_r(2:end-1);

save FORBILD_data img_70kev sino_mono_analytic sino_poly_analytic theta_deg numpix dx spectrum* energies* soft_vals bone_vals sino_mono_noisy sino_poly_noisy

%% generate so-called exact data. MIRT system matrix is used to generate the data, and attenuation exactly follows the model of pSART
load FORBILD_data

dims = size(sino_poly_analytic);

%read in attenuation coefficients
material_names = {'air','soft','bone'};
mass_atten=xray_read_atten(material_names,[1:energies(end)]);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
material_att = [air_vals'; soft_vals'; bone_vals'];

sino_mono_mirt = zeros(dims); sino_poly_mirt = sino_mono_mirt;
%create system matrix one view at a time (to prevent creating overly large matrix)
ig = image_geom('nx',numpix,'dx',dx);
for j = 1:numel(theta_deg)
    disp(j)
    sg = sino_geom('par','orbit_start',theta_deg(j),'orbit', 360, 'nb',dims(1),'na',1,'dr',dx,'strip_width','d');
    A = Gtomo2_strip(sg,ig);
    sino_mono_mirt(:,j) = exp(-A*img_70kev(:)); 
    sino_poly_mirt(:,j) = project_poly(img_70kev,Eref,energies,spectrum,material_att,material_names,A,[dims(1) 1],false);
end


save FORBILD_data img_70kev sino_mono_analytic sino_poly_analytic theta_deg numpix dx spectrum* energies* soft_vals bone_vals sino_mono_noisy sino_poly_noisy sino_poly_mirt sino_mono_mirt