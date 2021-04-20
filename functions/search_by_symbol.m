function [probe_index, probe_names] = search_by_symbol(T, symbol)
%SEARCH_BY_SYMBOL returns the probes positions of a gene symbol.
%   [PROBE_ID] = structure_id(T, NAMES) provided with the 
%   table of Probes, SEARCH_BY_SYMBOL matches its variable "gene_symbol" 
%   with the gene symbol SYMBOL. Each gene is probed with different
%   sequences, so there might me multiple probes for the gene of interest.
%   SEARCH_BY_SYMBOL returns the positions of the probes PROBE_INDEX, 
%   which can then beused to index the rows of the expression matrix. 
%   PROBE_NAMES contains the name of the corresponding probe. 
%   The function returns -1 if no gene matches the query.
mask = strcmpi(T.gene_symbol, symbol);
if ~any(mask)
    probe_index = -1;
    probe_names = [];
else
    probe_index = find(mask);
    probe_names = T.probe_name(mask);
end
