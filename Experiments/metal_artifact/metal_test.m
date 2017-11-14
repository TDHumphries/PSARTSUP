%% metal_test.m
% Code to generate reconstructions shown in "Superiorized polyenergetic reconstruction algorithm 
% for reduction of metal artifacts in CT images" by Humphries & Gibali (http://faculty.uwashington.edu/thumphri/HG17.pdf)
%
% Contact T. Humphries (thumphri@uw.edu) with questions

addpath(genpath('../../core'));
addpath('../../../Libraries/fessler/irt');

setup;
kVp = 130;
load(sprintf('metal_test_data_%d.mat',kVp));

win = [0.18 0.23];

%% reconstruction parameters

%geometric parameters
geom_params.g = 'par';
geom_params.arc = 180;
geom_params.start_angle = theta(1);
geom_params.numpix = numpix;
geom_params.dx = dx;

%reconstruction parameters (some to be filled in later)
recon_params.x_true = img_ref;
recon_params.x0 = zeros(numpix);
recon_params.beta = 1;
recon_params.use_weights = 0;
max_its = 32;   %store separately as the recon_params field will be overwritten for superiorized versions
recon_params.ns = 12;

%polyenergetic parameters
xrs=xray_read_spectra(sprintf('poly1,[%d]',kVp)); % use for I values
poly_params.energies = xrs.en;
poly_params.spectrum = xrs.sp/sum(xrs.sp);
poly_params.Eref = Eref;
poly_params.material_names = {'air','soft','bone','titanium'};

mass_atten=xray_read_atten(poly_params.material_names,poly_params.energies);
mass_dens = xray_read_dens(poly_params.material_names);

for j = 1:numel(poly_params.material_names)
    mass_atten(:,j) = mass_atten(:,j)*mass_dens(j);
end
poly_params.material_att = mass_atten';
poly_params.use_iodine = false;

%superiorization parameters
sup_params.gamma = 0.9995;
sup_params.alpha_init = 1;
sup_params.N = 40;
sup_params.objfun = @grad_TV;


%% Reconstructions at multiple count levels
counts = 1e5*[1 2 5 10];
for j = 1:numel(counts)
    fprintf('counts = %f\n',counts(j));
    sino_poly_noisy = 1e12/counts(j)*imnoise(counts(j)/1e12*sino_poly,'poisson'); sino_poly_noisy = max(sino_poly_noisy,1/counts(j));
    
    sino_poly_wc_noisy = soft_tissue_correction(sino_poly_noisy,energies,spectrum,soft_vals,Eref);
    
    recon_params.use_weights = 0; recon_params.max_its = max_its;
    fprintf('SART\n');[recon_mono1, ~, res1m, err1] =      pSART_sup(sino_poly_wc_noisy,geom_params,recon_params);
    sup_params.epsilon_target = res1m(end); recon_params.max_its = 1e6;
    fprintf('SART-sup\n');[recon_mono2, ~, res2m, err2] =  pSART_sup(sino_poly_wc_noisy,geom_params,recon_params,[],sup_params);
    recon_params.use_weights = 1; recon_params.max_its = max_its;
    fprintf('WSART\n');[recon_mono3, ~, res3m, err3] =      pSART_sup(sino_poly_wc_noisy,geom_params,recon_params);
    sup_params.epsilon_target = res3m(end);  recon_params.max_its = 1e6;
    fprintf('WSART-sup\n');[recon_mono4, ~, res4m, err4] = pSART_sup(sino_poly_wc_noisy,geom_params,recon_params,[],sup_params);

    recon_params.use_weights =0; recon_params.max_its = max_its;
    fprintf('pSART\n');[recon_poly1, ~, res1, err1] =      pSART_sup(sino_poly_noisy,geom_params,recon_params, poly_params);
    sup_params.epsilon_target = res1(end); recon_params.max_its = 1e6;
    fprintf('pSART-sup\n');[recon_poly2, ~, res2, err2] =  pSART_sup(sino_poly_noisy,geom_params,recon_params, poly_params,sup_params);
    recon_params.use_weights = 1; recon_params.max_its = max_its;
    fprintf('WpSART\n');[recon_poly3, ~, res3, err3] = pSART_sup(sino_poly_noisy,geom_params,recon_params, poly_params);
    sup_params.epsilon_target = res3(end); recon_params.max_its = 1e6;
    fprintf('WpSART-sup\n');[recon_poly4, ~, res4, err4] = pSART_sup(sino_poly_noisy,geom_params,recon_params, poly_params,sup_params);
    save(sprintf('results/recons_counts%f.mat',counts(j)),'recon_mono*','recon_poly*','res*','sino_poly_noisy');
end

%%
counts = 1e5*[1 2 5 10];
for j = 1:numel(counts)
    figure(j); colormap gray;
    load(sprintf('results/recons_counts%f.mat',counts(j)));
    subplot(2,4,1);imagesc(recon_mono1,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('SART');
    subplot(2,4,2);imagesc(recon_mono2,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('SART-TV');
    subplot(2,4,3);imagesc(recon_mono3,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WSART');
    subplot(2,4,4);imagesc(recon_mono4,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WSART-TV');
    subplot(2,4,5);imagesc(recon_poly1,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('pSART');
    subplot(2,4,6);imagesc(recon_poly2,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('pSART-TV');
    subplot(2,4,7);imagesc(recon_poly3,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WpSART');
    subplot(2,4,8);imagesc(recon_poly4,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WpSART-TV');
    figure(5); subplot(1,4,j); imagesc(-log(sino_poly_noisy)); axis image; colorbar; colormap gray;
    
    %print results
    fprintf('Count rate: %f\n',counts(j));
    TV = sup_params.objfun(recon_mono1); fprintf('SART:\t       Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res1m),res1m(end),TV);
    TV = sup_params.objfun(recon_mono2); fprintf('SART-TV:\t    Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res2m),res2m(end),TV);
    TV = sup_params.objfun(recon_mono3); fprintf('wSART:\t      Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res3m),res3m(end),TV);
    TV = sup_params.objfun(recon_mono4); fprintf('wSART-TV:\t   Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res4m),res4m(end),TV);
    
    TV = sup_params.objfun(recon_poly1); fprintf('pSART:\t       Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res1),res1(end),TV);
    TV = sup_params.objfun(recon_poly2); fprintf('pSART-TV:\t    Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res2),res2(end),TV);
    TV = sup_params.objfun(recon_poly3); fprintf('wpSART:\t      Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res3),res3(end),TV);
    TV = sup_params.objfun(recon_poly4); fprintf('wpSART-TV:\t   Its=%d \t epsilon = %4.2f \t TV = %5.1f\n',numel(res4),res4(end),TV);
  
end

%% NOISELESS TEST

    recon_params.use_weights =0; recon_params.max_its = max_its;
    fprintf('SART\n');[recon_mono1, ~, res1m, err1] =      pSART_sup(sino_poly_wc,geom_params,recon_params);
    sup_params.epsilon_target = res1m(end); recon_params.max_its = 1e6;
    fprintf('SART-sup\n');[recon_mono2, ~, res2m, err2] =  pSART_sup(sino_poly_wc,geom_params,recon_params,[],sup_params);
    recon_params.use_weights = 1; recon_params.max_its = max_its;
    fprintf('WSART\n');[recon_mono3, ~, res3m, err3] =      pSART_sup(sino_poly_wc,geom_params,recon_params);
    sup_params.epsilon_target = res3m(end); recon_params.max_its = 1e6;
    fprintf('WSART-sup\n');[recon_mono4, ~, res4m, err4] = pSART_sup(sino_poly_wc,geom_params,recon_params,[],sup_params);

    recon_params.use_weights =0; recon_params.max_its = max_its;
    fprintf('pSART\n');[recon_poly1, ~, res1, err1] =      pSART_sup(sino_poly,geom_params,recon_params, poly_params);
    sup_params.epsilon_target = res1(end); recon_params.max_its = 1e6;
    fprintf('pSART-sup\n');[recon_poly2, ~, res2, err2] =  pSART_sup(sino_poly,geom_params,recon_params, poly_params,sup_params);
    recon_params.use_weights = 1;  recon_params.max_its = max_its;
    fprintf('WpSART\n');[recon_poly3, ~, res3, err3] = pSART_sup(sino_poly,geom_params,recon_params, poly_params);
    sup_params.epsilon_target = res3(end); recon_params.max_its = 1e6;
    fprintf('WpSART-sup\n');[recon_poly4, ~, res4, err4] = pSART_sup(sino_poly,geom_params,recon_params, poly_params,sup_params);
    save('results/recons_noiseless.mat','recon_mono*','recon_poly*','res*','sino_poly_noisy');
    
    %%
    load('results/recons_noiseless.mat'),'recon_mono*','recon_poly*','res*','sino_poly_noisy');
    figure(11); colormap gray;
    subplot(2,4,1);imagesc(recon_mono1,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('SART');
    subplot(2,4,2);imagesc(recon_mono2,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('SART-TV');
    subplot(2,4,3);imagesc(recon_mono3,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WSART');
    subplot(2,4,4);imagesc(recon_mono4,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WSART-TV');
    subplot(2,4,5);imagesc(recon_poly1,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('pSART');
    subplot(2,4,6);imagesc(recon_poly2,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('pSART-TV');
    subplot(2,4,7);imagesc(recon_poly3,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WpSART');
    subplot(2,4,8);imagesc(recon_poly4,win); axis image; set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','bold'); title('WpSART-TV');