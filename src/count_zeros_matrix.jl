using Nemo

"""
    count_zero_rows(A)

Input:
- `A`: a matrix

The function counts the number of zero rows in `A`.

Output:
- `zero_rows`: the number of zero rows in `A`
"""
function count_zero_rows(A)
    return length([row for row in 1:n if iszero(A[row, :])])
end

"""
    count_zero_columns(A)

Input:
- `A`: a matrix

The function counts the number of zero columns in `A`.

Output:
- `zero_columns`: the number of zero columns in `A`
"""
function count_zero_columns(A)
    return length([col for col in 1:m if iszero(A[:, col])])
end
