function M = table_to_matrix(T)
    %  TABLE_TO_MATRIX Transform the table of expression in a matrix
    %  M = table_to_matrix(T) Takes the table T and returns the matrix M,
    %  without the first column (which contains probe indeces). 
M = T{:,2:end};