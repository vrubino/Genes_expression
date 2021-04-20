function mask = get_mask(all_ids, target_ids)
    %   get_masks return a mask to index the expression patterns.
    %   MASK = get_mask(ALL_IDS, TARGET_IDS) Provided with the list of all 
    %   ids in the 1xN vector ALL_IDS and the 1xM vector TARGET_IDS,
    %   returns a logical mask NxM mask; each row i is then 1xN vector
    %   which can be used to index the expression pattern and extract all
    %   the rows corresponding to the area i.
    arguments
        all_ids (:,1) {isnumeric}
        target_ids (:,1) {isnumeric}
    end
    mask = all_ids == target_ids';