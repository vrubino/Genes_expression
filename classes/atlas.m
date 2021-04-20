classdef atlas < handle
    properties
        Donors
        Probes
        IDs
        Names
        CommonStructuresId
        CommonStructuresNames
        CommonSignificantStructuresId
        CommonSignificantStructuresNames
        Paths
        PValues
    end
    
    methods
        
        % Constructor
        function obj = atlas(P, O)
            obj.Probes = P(:,[1 2 4 5]);
            obj.IDs = O.id;
            obj.Names = O.name;
            obj.Paths = O.structure_id_path;
        end
        
        % Add a donor
        function addDonor(obj, E, S)
            child = donor(E, S); % Construct the object for the donor
            obj.Donors = [obj.Donors child]; % Add the donor to the dataset
            commonStructures(obj, child.StructuresId);
        end
        
        function commonStructures(obj, ids)

            % Add common structures to the list
            if isempty(obj.CommonStructuresId)
                obj.CommonStructuresId = ids;
            else
                obj.CommonStructuresId = intersect(obj.CommonStructuresId,ids);
            end
            obj.CommonStructuresNames = id_to_name(obj.IDs, ...
                    obj.Names, obj.CommonStructuresId);
        end
        
        % Save the dataset for later use
        function store(obj, fileName)
            fh = @(x) inputname(1); % Get the name of any variable
            save(fileName, fh(obj));
        end
        
        % Slice the donor expressions by row (probe)
        function slice(obj, index)
            for ii = 1:numel(obj.Donors)
                slice(obj.Donors(ii), index)
            end
        end
        
        % Slice the expressions by columns (areas)
        function ROI(obj, id)
            if isnumeric(id)
                id = num2str(id);
            end
            
            obj.CommonStructuresId = []; % Reset the list of common structures
            paths = cell(numel(obj.Donors));
            for ii = 1:numel(obj.Donors)
                tmp_mask = obj.Donors(ii).NativeRowsId == obj.IDs';
                
                for jj = 1:size(tmp_mask,1)
                    paths{ii}(jj) = obj.Paths(tmp_mask(jj,:));
                end
                
                mask = ~cellfun(@(x) isempty(x), regexp(paths{ii}, id));
                obj.Donors(ii).NativeRowsId = obj.Donors(ii).NativeRowsId(mask);
                obj.Donors(ii).NewRowsId = obj.Donors(ii).NewRowsId(mask);
                obj.Donors(ii).Expr = obj.Donors(ii).Expr(:,mask);
                obj.Donors(ii).Mask = get_mask(obj.Donors(ii).NewRowsId, ...
                    unique(obj.Donors(ii).NewRowsId, 'stable')); 
                commonStructures(obj, obj.Donors(ii).NewRowsId);
            end
        end
        
        % Compute differential stability
        function DSs = DS(obj)
            n_probes = size(obj.Donors(1).Expr,1);
            
            % Collapse columns: as a single area is surveyed by multiple
            % samples, to compute DS we need to first average across
            % samples of the same areas
            avgExpr = cell(1,numel(obj.Donors));
            for ii = 1:numel(obj.Donors)
                for k = 1:size(obj.Donors(ii).Mask,2)
                    avgExpr{ii}(:,k) = mean(obj.Donors(ii).Expr(:,obj.Donors(ii).Mask(:,k)),2);
                end
            end
            
            % Compute NxN correlation matrix, where each (i,j) entry is the
            % correlation between donors i and j, across all areas
            [IA,IB] = deal(cell(1, numel(avgExpr)));
            rhos = cell(1,n_probes);
            for ii = 1:numel(obj.Donors)
                for jj = 1:numel(obj.Donors)
                    % Common areas
                    [~,IA{ii}{jj},IB{ii}{jj}] = intersect(obj.Donors(ii).StructuresId, ...
                        obj.Donors(jj).StructuresId);                  
                    
                    for k = 1:n_probes
                        rhos{k}(ii,jj) = corr(avgExpr{ii}(k,IA{ii}{jj})', ...
                            avgExpr{jj}(k,IB{ii}{jj})');
                    end
                end
            end
            
            % Retain only the upper triangle of the symmetric matrix
            % (discarding the meaningless diagonal)
            DSs = cellfun(@(x) mean(1-squareform(round(1-x,4))), rhos); 
            % First subtract from 1 to obtain distance matrix to use
            % squareform and then subtract again from 1
            [~,best_probe] = max(DSs);
            slice(obj, best_probe)
            
        end
        
        % Change Parcellation
        % Use the function handle FH to map each id onto another
        function changeParcellation(obj, fh)
            obj.CommonStructuresId = [];
            for ii = 1:numel(obj.Donors)
                NewIDs = fh(obj.IDs, obj.Paths, obj.Donors(ii).NativeRowsId);
                substituteStructures(obj.Donors(ii), NewIDs);
                commonStructures(obj, NewIDs);
            end
        end
        
        % Test
        function test(obj, correction, a)
            for ii = 1:numel(obj.Donors)
                % Perform the statistical test
                t_test(obj.Donors(ii), correction);
                
                % Update each donor with the list of significant structures
                whichStructures(obj.Donors(ii),a);
                
                % Get the list of common significant structures
                if ii == 1
                    CSS = obj.Donors(ii).SignificantStructures; % At the first 
                    % iteration, instatiate the variable with all
                    % significant structures of the first donor
                else
                    CSS = intersect(obj.Donors(ii).SignificantStructures, CSS);
                end
                obj.CommonSignificantStructuresId = CSS;
                obj.CommonSignificantStructuresNames = id_to_name(obj.IDs, ...
                    obj.Names, obj.CommonSignificantStructuresId);
            end
        end
        
        % Stack p-values
        function stack(obj)
            p = zeros(numel(obj.Donors), numel(obj.CommonStructuresId));
            for ii = 1:numel(obj.Donors)
                tmp_mask = obj.Donors(ii).StructuresId == obj.CommonStructuresId';
                for jj = 1:size(tmp_mask,2)
                    p(ii,jj) = obj.Donors(ii).CorrectedPValues(:,tmp_mask(:,jj));
                end  
            end
            obj.PValues = p;
        end
        
        function plot(obj, fileName, cmap)
            mask = cellfun(@(x) ~isempty(x), obj.CommonStructuresNames); 
            % Get indeces of structures with no names
            h = imagesc(obj.PValues(:,mask));
            xticks(1:numel(obj.CommonStructuresId(mask)))
            xtickangle(45)
            colormap(cmap)
            xticklabels(obj.CommonStructuresNames(mask))
            saveas(h, fileName)

        end
    end
end