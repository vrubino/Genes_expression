function names = id_to_name(all_ids, all_names, id)
    %   ID_TO_NAME retrieve the name of a structure by its id.
    %   NAMES = id_to_name(ALL_IDS, ALL_NAMES, ID): Provided with the list of ids
    %   ALL_IDS and of names ALL_NAMES, returns a cell NAMES containing the
    %   name of each of the ids contained in the vector ID of any
    %   dimensions.
    
names = cell(size(id));
for r = 1:size(id,1)
    for c = 1:size(id,2)
        if any(all_ids == id(r,c))
            names{r,c} = all_names{all_ids == id(r,c)};
        end
    end
end