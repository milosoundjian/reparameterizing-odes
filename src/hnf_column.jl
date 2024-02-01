using Nemo

"""
    hnf_column(A)

Input:
- `A`: a matrix of integers

The function calculates the Hermite column normal form of `A`.

Output:
- `H` the hermite column normal form of `A`
"""
# Gleb: Julia allows type annotations for function arguments, they typically
# make code more explanatory, here it should be A::ZZMatrix
function hnf_column(A)
    (m, n) = size(A)

    # make a copy of A
    A1 = deepcopy(A)

    # reverse the first m columns
    if n > m
        A1[:, 1:m] = A1[:, m:-1:1]
    else
        A1[:, 1:n] = A1[:, n:-1:1]
    end

    # reverse the rows
    A1 = A1[m:-1:1, :]

    H = transpose(hnf(transpose(A1)))

    # reverse the first m columns
    H[:, 1:m] = H[:, m:-1:1]

    # reverse the rows
    H = H[m:-1:1, :]

    return H
end

"""
    hnf_with_transform_column(A)

Input:
- `A`: a matrix of integers

The function calculates the Hermite column normal form of `A` and the transform matrix.

Output:
- `H`: the hermite column normal form of `A`
- `V`: The transform matrix. `A * V = H` and `V` is unimodular
"""
function hnf_with_transform_column(A)
    m = size(A)[1] # number of rows
    n = size(A)[2] # number of columns

    # make a copy of A
    A1 = deepcopy(A)
    # A1 = A

    # reverse the first m columns
    # if n > m
    #     A1[:, 1:m] = A1[:, m:-1:1]
    # else
    #     A1[:, 1:n] = A1[:, n:-1:1]
    # end

    # reverse the rows
    A1 = A1[m:-1:1, :]

    H, V = hnf_with_transform(transpose(A1))

    H = transpose(H)

    # reverse the first m columns
    if n > m
        H[:, 1:m] = H[:, m:-1:1]
    else
        H[:, 1:n] = H[:, n:-1:1]
    end
    # reverse the rows
    H = H[m:-1:1, :]

    # transpose the transform
    V = transpose(V)

    # reverse the first m columns
    if n > m
        V[:, 1:m] = V[:, m:-1:1]
    else
        V[:, 1:n] = V[:, n:-1:1]
    end

    return H, V
end
