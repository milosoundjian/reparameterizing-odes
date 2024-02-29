include("../src/hnf_column_normal.jl")

K = matrix(ZZ, [[1, 1], [0, -1], [1, 1], [0, 1]])

H_K, U_K = hnf_with_transform(K)

H_K

U_K

r = count_zero_rows(H_K)

A = U_K[end-r+1:end, :]

_, V = hnf_with_normal_transform_column(A)

V

W = inv(V)