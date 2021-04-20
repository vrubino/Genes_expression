function ancestors = get_ancestor(all_ids, paths, K, target_ids)
    %   GET_ANCESTOR retrieve the ancestors of the input areas.
    %   [ANCESTORS] = get_ancestor(ALL_IDS, PATHS, K, TARGET_IDS)
    %   Provided with the list of ids ALL_IDS and their respective paths in
    %   PATHS, it returns the ANCESTORS at the K-degree for the areas whose
    %   ids are in TARGET_IDS
    struct_paths = cellfun(@(x) regexp(x, '\d{4}', 'match'), ...
        paths, 'un', 0);
    K_ancestors = cellfun(@(x) str2num(x{min(K, numel(x))}), struct_paths);
    
    ancestors = zeros(numel(target_ids),1);
    for ii = 1:numel(target_ids)
        ancestors(ii) = K_ancestors(target_ids(ii) == all_ids);
    end
    
    
