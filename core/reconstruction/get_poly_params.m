function [energies,spectrum,Eref, material_att, material_names, mask, use_iodine] = get_poly_params(params)
    
    fprintf('Reading polyenergetic parameters...\n');
    
    %energies
    assert(isfield(params,'energies'),'Must specify poly_params.energies for polyenergetic reconstruction')
    energies = params.energies(:);
    ne = numel(energies);
  
    %spectrum
    assert(isfield(params,'spectrum'),'Must specify poly_params.spectrum for polyenergetic reconstruction')
    spectrum = params.spectrum(:);
    assert(numel(spectrum)==ne,'size of energies and spectrum are not the same');
    
    %reference energy
    if ~isfield(params,'Eref')
        fprintf('\tReference energy not specified in poly_params.Eref -- using 70 keV\n');
       	Eref = 70;
    else
        Eref = params.Eref;
    end    

    %material attenuation
    assert(isfield(params,'material_att'),'Must specify poly_params.material_att (M x E matrix) for polyenergetic reconstruction');

    material_att = params.material_att;
    assert(size(material_att,2)==ne,'Number of columns in poly_params.material_att must match size of poly_params.energies');
    nm = size(material_att,1);
    
    %using iodine
    if ~isfield(params,'use_iodine')
        fprintf('\tpoly_params.use_iodine not specified -- assuming no iodine\n');
       	use_iodine = false;
    else
        use_iodine = params.use_iodine;
    end
    
    %need material names and mask
    if use_iodine
        assert(isfield(params,'material_names'),'Must include material names in poly_params.material_names if using iodine');
        material_names = params.material_names;
        assert(numel(material_names) == nm,'Number of materials in material_names does not match material_att');
        
        assert(isfield(params,'mask'),'Must include binary mask in poly_params.mask if using iodine');
        mask = params.mask;
    else
        mask = []; material_names = {};
    end
end