using Nemo
include("hnf_column.jl")


A = matrix(ZZ, [[3, 3, 1, 4], [1, 2, 1, 3], [1, 1, 1, 2], [5, 3, 4, 2], [1, 1, 2, 3]])

H, V = hnf_with_transform(A)
V * A == H

hnf_column(A)

hnf_with_transform_column(hnf_with_transform_column(A)[1])
