function [geom, arc, start_angle, numpix, dx, fb_params] = get_geom_params(params)
    fprintf('Reading geometry parameters...\n');
   
    %check geometry and fan-beam parameters
    if ~isfield(params,'g')
        fprintf('\tNo geometry specified in geom_params.g -- assuming parallel\n');
        geom = 'par';
    else
        switch params.g                 %parallel beam
            case 'par'
                geom = params.g;
                fb_params = [];
            case 'fan'                  %fan beam -- need additional parameters
                geom = params.g;
                
                assert(isfield(params,'dso'),'Must specify source object distance (geom_params.dso) for fan beam')
                fb_params.dso = params.dso;
                
                assert(isfield(params,'dod'),'Must specify object detector distance (geom_params.dod) for fan beam')
                fb_params.dod = params.dod;
                
                assert(isfield(params,'ds'),'Must specify detector spacing (geom_params.ds) for fan beam')
                fb_params.ds = params.ds;
                
            otherwise                   %invalid input
                error('Unrecognized option: %s in params.g. Valid options are ''par'' or ''fan''.',params.g)
        end
    end
    
    %check other parameters
    if ~isfield(params,'arc')                  %acquisition arc
        fprintf('\tArc not specified in geom_params.arc -- assuming 360 degrees\n');
        arc = 360;
    else
        arc = params.arc;
    end
    
    if ~isfield(params,'start_angle')          %start angle
        fprintf('\tStart angle not specified in geom_params.start_angle -- assuming 0 degrees\n');
        start_angle =0;
    else
        start_angle = params.start_angle;
    end
    
    assert(isfield(params,'numpix'),'Must specify number of pixels (geom_params.numpix)');
    numpix = params.numpix;
    
    assert(isfield(params,'dx'),'Must specify pixel width (geom_params.dx)');
    dx = params.dx;

    
end