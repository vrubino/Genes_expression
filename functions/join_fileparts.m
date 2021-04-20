function file_path = join_fileparts(varargin)
%JOIN_FILEPARTS Joins into char vectors file parts, using / as delimiter.
%   NEWSTR = join_fileparts(S1,S2,...) joins S1,S2,... with / as delimiter.
%
file_path = strjoin(varargin(cellfun(@(x) ~isempty(x), varargin)),'/');
end