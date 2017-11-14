%% exp1_sparseview.m
% Code used for sparse-view experiments in "Superiorized algorithm for reconstruction
% of CT images from sparse-view and limited-angle polyenergetic data" by
% Humphries, Winn and Faridani. https://arxiv.org/abs/1701.03396
%
% Contact T. Humphries (thumphri@uw.edu) with questions

%%
addpath(genpath('../../core'));
addpath('../../../Libraries/fessler/irt');  %put your own path to IRT here

setup;

%%
win = [0.195 0.215];                %display window for colormap (mu values)
da = [1 2 4 5 10 20];               %spacing between views. da = 1 corresponds to complete data, da = 2 skips every other view, etc.
numits = [18 36 72 90 180 360];    %number of iterations to run pSART for consistent data
ns = [120 60 30 24 12 6];           %number of subsets of projection data

load FORBILD_data.mat


%% parameters for reconstruction

%geometric parameters
geom_params.g = 'par';
geom_params.arc = 180;
geom_params.start_angle = theta(1)+90;   %theta is in radon transform convention, need to convert to angular (view) convention
geom_params.numpix = numpix;
geom_params.dx = dx;

%reconstruction parameters (some to be filled in later)
recon_params.x_true = img_ref;
recon_params.x0 = zeros(numpix);
recon_params.beta = 1;
recon_params.use_weights = 0; %not used in these experiments

%polyenergetic parameters
poly_params.energies = energies;
poly_params.spectrum = spectrum;
poly_params.Eref = Eref;
poly_params.material_names = {'air','soft','bone'};
mass_atten=xray_read_atten(poly_params.material_names,energies);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];
poly_params.use_iodine = false;

%superiorization parameters
sup_params.gamma = 0.999;
sup_params.alpha_init = 1;
sup_params.N = 20;
sup_params.objfun = @grad_TV;


%% Polyenergetic experiments (exact) -- Top row of Figure 4, left section of Table 1
for j = 1:6
    fprintf('%d views\n',1440/da(j));
    data = sino_poly_mirt(:,1:da(j):1440); theta = theta_deg(1:da(j):1440);
    
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
   
    sup_params.epsilon_target = res_pSART(end);
    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    recon_params.max_its = 1e6; %run until epsilon_target is reached
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%dviews_poly_exact.mat',1440/da(j)),'recon*','res_*','err_*');
end

%% Polyenergetic experiments (inconsistent data, count rate of 4e6) -- Second row of figure 4

numits = [9 18 36 45 90 180];    %reduce number of iterations for inconsistent data

poly_params.energies = energies_r;
poly_params.spectrum = spectrum_r;
mass_atten=xray_read_atten(poly_params.material_names,energies_r);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];


for j = 1:6
    fprintf('%d views\n',1440/da(j));
    data = sino_poly_noisy(:,1:da(j):1440); theta = theta_deg(1:da(j):1440);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
   
    sup_params.epsilon_target = res_pSART(end);
    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    recon_params.max_its = 1e6; %run until epsilon_target is reached
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%dviews_poly_inconsistent.mat',1440/da(j)),'recon*','res_*','err_*');
end

%% Polyenergetic experiments (inconsistent data, count rate of 1e6)  -- Third row of figure 4

counts = 1e6;
sino_mono_noisy = 1e12/counts*imnoise(counts/1e12*sino_mono_analytic,'poisson');
sino_poly_noisy = 1e12/counts*imnoise(counts/1e12*sino_poly_analytic,'poisson');


for j = 1:6
    fprintf('%d views\n',1440/da(j));
    data = sino_poly_noisy(:,1:da(j):1440); theta = theta_deg(1:da(j):1440);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
   
    sup_params.epsilon_target = res_pSART(end);
    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    recon_params.max_its = 1e6; %run until epsilon_target is reached
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%dviews_poly_inconsistent_%fcounts.mat',1440/da(j),counts),'recon*','res_*','err_*');
end

%% Polyenergetic experiments (inconsistent data, 80kVp spectrum)  -- fourth row of figure 4
load FORBILD_data_80kVp.mat

recon_params.x_true = img_50kev;
poly_params.energies = energies_r;
poly_params.spectrum = spectrum_r;
poly_params.Eref = 50;

mass_atten=xray_read_atten(poly_params.material_names,energies_r);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];


for j = 1:6
    fprintf('%d views\n',1440/da(j));
    data = sino_poly_noisy(:,1:da(j):1440); theta = theta_deg(1:da(j):1440);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
   
    sup_params.epsilon_target = res_pSART(end);
    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    recon_params.max_its = 1e6; %run until epsilon_target is reached
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%dviews_poly_inconsistent_80kVp_test.mat',1440/da(j)),'recon*','res_*','err_*');
end

%% generate figure
xmin = 130; xmax = 670; xsize = xmax - xmin;
ymin = 70; ymax = 730; ysize = ymax - ymin; 
F = zeros(ysize*3,xsize*6); F2 = zeros(ysize,xsize*6);
nviews = [1440 288 144];
for j = 1:3
    load(sprintf('results/FORBILD_%dviews_poly_exact.mat',nviews(j)));
    F(1:ysize, ((2*j-2)*xsize + 1) : (2*j-1)*xsize) = recon_pSART(ymin+1:ymax,xmin+1:xmax);
    F(1:ysize, ((2*j-1)*xsize + 1) : 2*j*xsize) = recon_sup(ymin+1:ymax,xmin+1:xmax);
    
    load(sprintf('results/FORBILD_%dviews_poly_inconsistent.mat',nviews(j)));
    F(ysize+1:2*ysize, ((2*j-2)*xsize + 1) : (2*j-1)*xsize) = recon_pSART(ymin+1:ymax,xmin+1:xmax);
    F(ysize+1:2*ysize, ((2*j-1)*xsize + 1) : 2*j*xsize) = recon_sup(ymin+1:ymax,xmin+1:xmax);
    
    load(sprintf('results/FORBILD_%dviews_poly_inconsistent_%fcounts.mat',nviews(j),1e6));
    F(2*ysize+1:3*ysize, ((2*j-2)*xsize + 1) : (2*j-1)*xsize) = recon_pSART(ymin+1:ymax,xmin+1:xmax);
    F(2*ysize+1:3*ysize, ((2*j-1)*xsize + 1) : 2*j*xsize) = recon_sup(ymin+1:ymax,xmin+1:xmax);
    
    load(sprintf('results/FORBILD_%dviews_poly_inconsistent_80kVp.mat',nviews(j)));
    F2(1:ysize, ((2*j-2)*xsize + 1) : (2*j-1)*xsize) = recon_pSART(ymin+1:ymax,xmin+1:xmax);
    F2(1:ysize, ((2*j-1)*xsize + 1) : 2*j*xsize) = recon_sup(ymin+1:ymax,xmin+1:xmax);
end
figure(22); imagesc(F,[0.195 0.215]); axis image; colormap gray;
set(gca,'FontSize',16,'FontWeight','bold');
set(gca,'XTick',xsize*(0:6) + xsize/2,'YTick',ysize/2*[1 3 5],'XTickLabel',{'pSART','pSART-TV','pSART','pSART-TV','pSART','pSART-TV'},'YTickLabel',{'C','I-1','I-2'},'XAxisLocation','top');
figure(33); imagesc(F2,[0.235 0.255]); axis image; colormap gray;
set(gca,'XTick',[],'YTick',ysize/2,'YTickLabel','I-3');
