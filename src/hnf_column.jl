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
    r = rank(A)

    # make a copy of A
    A1 = deepcopy(A)

    # reverse the rows
    A1 = A1[m:-1:1, :]

    H = transpose(hnf(transpose(A1)))

    # reverse the rows
    H = H[m:-1:1, :]
    
    # reverse the first r columns
    H[:, 1:r] = H[:, r:-1:1]

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
    (m, n) = size(A)
    r = rank(A)

    # make a copy of A
    A1 = deepcopy(A)

    # reverse the rows
    A1 = A1[m:-1:1, :]

    H, V = hnf_with_transform(transpose(A1))

    H = transpose(H)

    # reverse the rows
    H = H[m:-1:1, :]

    # reverse the first r columns
    H[:, 1:r] = H[:, r:-1:1]

    # transpose the transform
    V = transpose(V)

    # reverse the first r columns
    V[:, 1:r] = V[:, r:-1:1]

    return H, V
end
