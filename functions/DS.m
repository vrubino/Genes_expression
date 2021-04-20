function [DSs,rhos] = DS(Expr, IDs, mask)

n_probes = size(Expr{1},1);

%% 1. Pool together
aExpr = cellfun(@(x,y) collapse_columns(x,y), Expr, mask, 'un', 0);

% Pre-allocate
[IA,IB] = deal(cell(1, numel(aExpr)));
[IA{:}] = deal(cell(1, numel(aExpr)));
[IB{:}] = deal(cell(1, numel(aExpr)));
rhos = cell(1, n_probes);

% Loop across donors
for ii = 1:numel(IDs)
    for jj = 1:numel(IDs)
        
        % Get indexes of common structures between donors i and j
        [~,IA{ii}{jj},IB{ii}{jj}] = intersect(IDs{ii}, IDs{jj});
        
        for k = 1:n_probes
            rhos{k}(ii,jj) = corr(aExpr{ii}(k,IA{ii}{jj})', ...
                aExpr{jj}(k,IB{ii}{jj})');
        end
    end
end

DSs = cellfun(@(x) mean(squareform(round(x,4))), rhos);