using Nemo

# Verhulst Model of Logistic Growth
# dn/dt = r * n * (1 - n/k)

# Let us get the Laurent Polynomial
R, (r, k, t, n) = polynomial_ring(QQ, ["r", "k", "t", "n"])

F = [r * t - r * (1//k) * t * n]
F = [t * r * (1 - n // k)]

S, (t1, n1) = laurent_polynomial_ring(QQ, ["t1", "n1"])

result = evaluate(F[1], [one(S), one(S), t1, n1])

A = matrix(ZZ, [[1, 2, 3], [4, 5, 6], [7, 8, 9]])

H, U = hnf_with_transform(A)

H
U * A == H

A = matrix(ZZ, [[3, 3, 1, 4], [1, 2, 1, 3], [1, 1, 1, 2]])

hnf(A)

include("../src/hnf_column.jl")

hnf_column(A)

hnf_column(transpose(A))

B = matrix(ZZ, [[1, 1, 1, 1, 1], [2, 2, 2, 2, 2], [1, 2, 3, 4 ,5], [2, 4, 6, 8 ,10], [3, 6, 9, 12, 15]])

hnf(B)

hnf_column(B)

C = matrix(ZZ, [[1, 2, 3, 4, 5], [2, 4, 6, 8, 10], [3, 6, 9, 12, 15], [4, 8, 12, 16, 20], [5, 10, 15, 20, 25]])

hnf(C)

hnf_column(C)