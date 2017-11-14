function [] = print_params(use_poly,use_iodine,use_sup,use_weights) 

fprintf('\nBeginning pSART iteration\n\t');
    if (use_poly)
        if (use_iodine)
            fprintf('Using polyenergetic projection with iodine\n\t');
        else
            fprintf('Using polyenergetic projection without iodine\n\t');
        end
    else
        fprintf('Using monoenergetic projection\n\t');
    end
    
    if (use_sup)
        fprintf('Using superiorization\n\t')
    else
        fprintf('Not using superiorization\n\t');
    end
    
    if (use_weights)
        fprintf('Using weighting\n')
    else
        fprintf('Not using weighting\n');
    end
end