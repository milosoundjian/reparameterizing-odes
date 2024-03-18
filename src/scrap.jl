using Nemo

# Verhulst Model of Logistic Growth
# dn/dt = r * n * (1 - n/k)

# Let us get the Laurent Polynomial
R, (r, k, t, n) = polynomial_ring(QQ, ["r", "k", "t", "n"])

F = [r * t - r * (1 // k) * t * n]
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

B = matrix(ZZ, [[1, 1, 1, 1, 1], [2, 2, 2, 2, 2], [1, 2, 3, 4, 5], [2, 4, 6, 8, 10], [3, 6, 9, 12, 15]])

hnf(B)

hnf_column(B)

hnf(C)

hnf_column(C)
+

C = matrix(ZZ, [[1, 2, 3, 4, 0], [2, 4, 6, 8, 0], [3, 6, 9, 12, 0], [4, 8, 12, 16, 0], [0, 0, 0, 0, 0]])

H, V = hnf_with_transform_column(C)

A = matrix(ZZ, [[1, 2, 3], [4, 5, 6], [5, 7, 9]])

function create_random_integer_matrix(m, n)
    return matrix(ZZ, rand(1:100, m, n))
end

function fake_hnf_column(A)
    (m, n) = size(A)
    r = rank(A)

    # make a copy of A
    A1 = deepcopy(A)

    # reverse the first m columns
    # if n > m
    #     A1[:, 1:m] = A1[:, m:-1:1]
    # else
    #     A1[:, 1:n] = A1[:, n:-1:1]
    # end

    # reverse the rows
    # A1 = A1[m:-1:1, :]

    H = transpose(hnf(transpose(A1)))

    # reverse the first r columns
    # H[:, 1:r] = H[:, r:-1:1]

    # reverse the rows
    # H = H[m:-1:1, :]

    return H
end

A = matrix(ZZ, [[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])
H_true = matrix(ZZ, [[3, 2, 0, 0, 0], [0, 1, 0, 0, 0]])
H1 = hnf_column(A)
H2, V = hnf_with_transform_column(A)
V
A * V

include("hnf_column.jl")

A = matrix(ZZ, [[6, 0, -4, 1, 3], [4, 3, 1, -4, 3]])

hnf(A)

hnf_column(A)

hnf_with_transform_column(A)[2]

odepp = @ODEmodel(
    n'(t) = n * (r * (1 - n/K - k * (p / (n+d)))),
    p'(t) = s * p * (1 - h * (p / n))
)