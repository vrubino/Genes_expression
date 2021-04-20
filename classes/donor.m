classdef donor < handle
    properties
        Expr
        StructuresId
        NativeRowsId
        NewRowsId
        Mask
        CorrectedPValues
        SignificantStructures
    end
    
    methods
        % Constructor
        function obj = donor(E, S)
            obj.Expr = E;
            obj.NativeRowsId = S.structure_id;
            obj.NewRowsId = obj.NativeRowsId;
            obj.StructuresId = unique(obj.NewRowsId, 'stable');
            obj.Mask = get_mask(obj.NativeRowsId,obj.StructuresId);  
        end
        
        % Slice the matrix by rows (probes)
        function slice(obj, index)
            obj.Expr = obj.Expr(index,:);
        end
        
        % Change parcellation
        function substituteStructures(obj, ids)
            obj.NewRowsId = ids;
            obj.StructuresId = unique(obj.NewRowsId, 'stable');
            obj.Mask = get_mask(obj.NewRowsId,obj.StructuresId); 
        end
        
        % Perform t-test
        function t_test(obj, correction)
            p = zeros(size(obj.Expr,1),size(obj.Mask,2));
            cp = zeros(size(obj.Expr,1),size(obj.Mask,2));
            for k = 1:size(obj.Mask,2)
                for ii = 1:size(obj.Expr,1)
                    tmp_mask = obj.Mask(:,k);
                    [~,p(ii,k)] = ttest2(obj.Expr(ii,tmp_mask), obj.Expr(ii,~tmp_mask), ...
                                'tail', 'right', 'Vartype','unequal');
                end
            end
            
            for ii = 1:size(obj.Expr,1)
                cp(ii,:) = correction(p(ii,:));
            end
            obj.CorrectedPValues = cp;
        end
        
        % Which structures have significatly higher expression?
        function whichStructures(obj, a)
            ids = cell(size(obj.CorrectedPValues,1),1);
            for ii = 1:size(obj.CorrectedPValues,1)
                tmp_mask = obj.CorrectedPValues(ii,:) < a;
                ids{ii} = obj.StructuresId(tmp_mask);
            end

            obj.SignificantStructures = cell2mat(ids);
        end
        
        

    end
end