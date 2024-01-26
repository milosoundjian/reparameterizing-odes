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
    n = size(A, 1)
    zero_rows = 0
    for row in 1:n
        if all(A[row, :] .== 0)
            zero_rows += 1
        end
    end
    return zero_rows
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
    m = size(A, 2)
    zero_columns = 0
    for column in 1:m
        if all(A[:, column] .== 0)
            zero_columns += 1
        end
    end
    return zero_columns
end
