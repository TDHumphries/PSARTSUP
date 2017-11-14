%% exp1_limitedangle.m
% Code used for limited angle-view experiments in "Superiorized algorithm for reconstruction
% of CT images from sparse-view and limited-angle polyenergetic data" by
% Humphries, Winn and Faridani. https://arxiv.org/abs/1701.03396
%
% Contact T. Humphries (thumphri@uw.edu) with questions

addpath(genpath('../../core'));
addpath('../../../Libraries/fessler/irt');  %put your own path to IRT here

setup;

%%
load FORBILD_data.mat
win     = [0.195 0.215];                                        %display window for colormap (mu values)
ns      = [110 110 100 100 90 90 80 80];                        %number of subsets of projection data
numits  = [75 75 85 85 95 95 105 105];                          %number of iterations to run pSART for consistent data
starts  = [97.5 7.49 104.9 14.9 112.49 22.49 119.99 29.99];     % starting angle for each acquisition arc
extents = [165 165 150 150 135 135 120 120];                    %angular extent of each acquisition arc

material_list = {'air','soft','bone'};

% Original theta's range from -180 to 180. Change to range from 0 to 360.
% for easier calculation of ATV weights. Also sort sinograms to be in order
% from 0 to 360.
ind = theta_deg < 0;
theta_deg(ind) = theta_deg(ind)+360;
[theta_deg,ind] = sort(theta_deg);
sino_poly_mirt = sino_poly_mirt(:,ind);
sino_poly_analytic = sino_poly_analytic(:,ind);
sino_poly_noisy = sino_poly_noisy(:,ind);

%for j = 1:numel(theta_deg)
%    element = theta_deg(j);
%    if element < 0
%        theta_deg(j) = element + 360;
%    end
%end

%geometric parameters (some to be filled in later)
geom_params.g = 'par';
geom_params.numpix = numpix;
geom_params.dx = dx;

%reconstruction parameters (some to be filled in later)
recon_params.x_true = img_70kev;
recon_params.x0 = zeros(numpix);
recon_params.beta = 1;
recon_params.use_weights = 0;

%polyenergetic parameters
poly_params.energies = energies;
poly_params.spectrum = spectrum;
poly_params.Eref = Eref;
poly_params.material_names = {'air','soft','bone'};
mass_atten=xray_read_atten(poly_params.material_names,energies);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];
poly_params.use_iodine = false;


%parameters for pSART-sup
alphas  = [0 45 90 135];    %directions used for computation of ATV

%superiorization parameters
sup_params.gamma = 0.9999;
sup_params.alpha_init = 0.4;
sup_params.N = 100;


%% Polyenergetic (exact). Top row of Figure 6. only do 165 and 150 degree extents
for j = 1:4
    tic
    % Extract correct range of angles
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    
    ind = theta_deg >= start_extent & theta_deg < end_extent;
    data = sino_poly_mirt(:,ind); theta = theta_deg(ind);

    geom_params.arc = extents(j); geom_params.start_angle = theta(1);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
     
    fprintf('start: %.2f \t end: %.2f \n', start_extent, end_extent);
    fprintf('views: %.2f \t views divided by ns: %.2f\n', numel(theta), numel(theta)/ns(j));
    
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
    %load(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_exact.mat', start_extent, end_extent));
    
    sup_params.epsilon_target = res_pSART(end);
    %objective function depends on acquisition arc
    sup_params.objfun = @(x) grad_ATV(x,theta,alphas);
    recon_params.max_its = 1e6;

    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_exact.mat', start_extent, end_extent), 'recon*','res_*','err_*', 'theta', 'alphas');
    toc
end

%% Polyenergetic (inconsistent). Right part of Table 2 and top row of Figure 6
numits  = [15 15 17 17 19 19 21 21]; %fewer iterations for inconsistent data

poly_params.energies = energies_r;
poly_params.spectrum = spectrum_r;
mass_atten=xray_read_atten(poly_params.material_names,energies_r);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];

for j = 1:4
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    
    ind = theta_deg >= start_extent & theta_deg < end_extent;
    data = sino_poly_noisy(:,ind); theta = theta_deg(ind);
    
    geom_params.arc = extents(j); geom_params.start_angle = theta(1);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    
    fprintf('start: %.2f \t end: %.2f \n', start_extent, end_extent);
    fprintf('views: %.2f \t views divided by ns: %.2f\n', numel(theta), numel(theta)/ns(j));
    
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
    
    sup_params.epsilon_target = res_pSART(end);
    %objective function depends on acquisition arc
    sup_params.objfun = @(x) grad_ATV(x,theta,alphas);
    recon_params.max_its = 1e6;

    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent.mat', start_extent, end_extent), 'recon*','res_*','err_*', 'theta', 'alphas');
    
end

%% Polyenergetic experiments (inconsistent data, higher noise level) 

counts = 1e6;
sino_mono_noisy = 1e12/counts*imnoise(counts/1e12*sino_mono_analytic,'poisson');
sino_poly_noisy = 1e12/counts*imnoise(counts/1e12*sino_poly_analytic,'poisson');

for j = 1:4
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    
    ind = theta_deg >= start_extent & theta_deg < end_extent;
    data = sino_poly_noisy(:,ind); theta = theta_deg(ind);
    
    geom_params.arc = extents(j); geom_params.start_angle = theta(1);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    
    fprintf('start: %.2f \t end: %.2f \n', start_extent, end_extent);
    fprintf('views: %.2f \t views divided by ns: %.2f\n', numel(theta), numel(theta)/ns(j));
    
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
    
    sup_params.epsilon_target = res_pSART(end);
    %objective function depends on acquisition arc
    sup_params.objfun = @(x) grad_ATV(x,theta,alphas);
    recon_params.max_its = 1e6;

    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent_counts%f.mat', start_extent, end_extent,counts), 'recon*','res_*','err_*', 'theta', 'alphas');
    
end

%% Polyenergetic experiments (low kVp spectrum) 

load FORBILD_data_80kVp.mat

ind = theta_deg < 0;
theta_deg(ind) = theta_deg(ind)+360;
[theta_deg,ind] = sort(theta_deg);
sino_poly_noisy = sino_poly_noisy(:,ind);

recon_params.x_true = img_50kev;
poly_params.energies = energies_r;
poly_params.spectrum = spectrum_r;
poly_params.Eref = 50;

mass_atten=xray_read_atten(poly_params.material_names,energies_r);
air_vals = mass_atten(:,1)*1.2e-3; soft_vals = mass_atten(:,2)*1.06; bone_vals = mass_atten(:,3)*1.92; %multiply by densities to get LAC
poly_params.material_att = [air_vals'; soft_vals'; bone_vals'];

numits = [23 23 25 25];
Eref = 50;
for j = 1:4
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    
    ind = theta_deg >= start_extent & theta_deg < end_extent;
    data = sino_poly_noisy(:,ind); theta = theta_deg(ind);
    
    geom_params.arc = extents(j); geom_params.start_angle = theta(1);
    recon_params.max_its = numits(j); recon_params.ns = ns(j);
    
    fprintf('start: %.2f \t end: %.2f \n', start_extent, end_extent);
    fprintf('views: %.2f \t views divided by ns: %.2f\n', numel(theta), numel(theta)/ns(j));
    
    [recon_pSART,~,res_pSART,err_pSART] = pSART_sup_new(data,geom_params,recon_params,poly_params);
    
    sup_params.epsilon_target = res_pSART(end);
    %objective function depends on acquisition arc
    sup_params.objfun = @(x) grad_ATV(x,theta,alphas);
    recon_params.max_its = 1e6;

    fprintf('Calling pSART-sup; target %3.2f\n',sup_params.epsilon_target);
    [recon_sup,~,res_sup,err_sup] = pSART_sup_new(data,geom_params,recon_params,poly_params,sup_params);
    save(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent_80kVp.mat', start_extent, end_extent), 'recon*','res_*','err_*', 'theta', 'alphas');
end

%% ATV Values of True Phantom at various extents
fprintf('True ATV: %.4g \n', grad_ATV(img_70kev, theta_deg, alphas));
ATV_values = zeros(1, 4);
for j = 1:4
    tic
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    ind = theta_deg >= start_extent & theta_deg < end_extent;
    theta = theta_deg(ind);
    
    TV = grad_ATV(img_50kev, theta, alphas);
    ATV_values(j) = TV;
    fprintf('Limited Extent ATV: %.4g \n', TV);
end
fprintf('Avg ATV: %.4g \n', mean(ATV_values));


%% Display results
for j = 1:1
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    load(sprintf('../../experiments_new/PMBPaper/results/FORBILD_%.2fstart_%.2fend_poly_inconsistent.mat', start_extent, end_extent));
    fprintf('start: %.2f \t end: %.2f \t views: %.2f\n', start_extent, end_extent, numel(theta));
    TV1 = grad_ATV(recon_pSART, theta, alphas); TV2 = grad_ATV(recon_sup, theta, alphas);
    figure(j); imagesc(recon_pSART,win); axis image; colormap gray;
    figure(j*10); imagesc(recon_sup,win); axis image; colormap gray;
    title(['Start: ' num2str(start_extent)  '   end: ' num2str(end_extent)])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    fprintf('pSART: numits = %d eps = %.3g TV = %.4g \npSART-sup: numits = %d eps = %.3g TV = %.4g \n',numel(res_pSART),res_pSART(end),TV1,numel(res_sup),res_sup(end),TV2);
    % pause();
end

%% generate figure
F = zeros(800*3,800*8); F2 = zeros(800,800*8);
for j = 1:4
    start_extent = starts(j);
    end_extent = start_extent + extents(j);
    
    load(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_exact.mat', start_extent, end_extent));
    F(1:800, ((2*j-2)*800 + 1) : (2*j-1)*800) = recon_pSART;
    F(1:800, ((2*j-1)*800 + 1) : 2*j*800) = recon_sup;
  
    load(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent.mat', start_extent, end_extent));
    F(801:1600, ((2*j-2)*800 + 1) : (2*j-1)*800) = recon_pSART;
    F(801:1600, ((2*j-1)*800 + 1) : 2*j*800) = recon_sup;
  
    load(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent_counts%f.mat', start_extent, end_extent,1e6));
    F(1601:2400, ((2*j-2)*800 + 1) : (2*j-1)*800) = recon_pSART;
    F(1601:2400, ((2*j-1)*800 + 1) : 2*j*800) = recon_sup;    
    
    load(sprintf('results/FORBILD_%.2fstart_%.2fend_poly_inconsistent_80kVp.mat',start_extent, end_extent));
    F2(1:800, ((2*j-2)*800 + 1) : (2*j-1)*800) = recon_pSART;
    F2(1:800, ((2*j-1)*800 + 1) : 2*j*800) = recon_sup;
end
figure(22); imagesc(F,[0.195 0.215]); axis image; colormap gray;
set(gca,'XTick',[],'YTick',[]);
figure(33); imagesc(F2,[0.235 0.255]); axis image; colormap gray;
set(gca,'XTick',[],'YTick',[]);