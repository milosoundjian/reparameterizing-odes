include("hnf_column.jl")
include("count_zeros_matrix.jl")

function hnf_with_normal_transform_column(A)
    # if A is tall, then we just perform hnf_with_transform_column
    if size(A)[1] > size(A)[2]
        return hnf_with_transform_column(A)
    end

    H = hnf_column(A)

    r = size(A)[1] # number of rows
    n = size(A)[2] # number of columns

    # if A is fat, then we stack a square matrix on top of a
    B = vcat(identity_matrix(ZZ, n), A)

    # We then get the HNF of this matrix
    H_star, V_star = hnf_with_transform_column(B)

    # V_i is the last r columns of H_star with the bottom r rows removed
    V_i = H_star[1:n, n-r+1:n]

    # V_n is the first n-r columns of H_star with the bottom r rows removed
    V_n = H_star[1:n, 1:n-r]

    V = hcat(V_i, V_n)

    return H, V
end
