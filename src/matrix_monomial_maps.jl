using Nemo

function split_matrix_pos_and_neg(A)
    # Splits a matrix A into two parts A_pos and A_neg
    # where A_pos contains the positive entries of A and 0s elsewhere
    # and A_neg contains the negative entries of A and 0s elsewhere
    n, m = size(A)
    A_pos = deepcopy(A)
    A_neg = deepcopy(A)

    for i in 1:n
        for j in 1:m
            if A[i, j] < 0
                A_pos[i, j] = 0
            else
                A_neg[i, j] = 0
            end
        end
    end

    return A_pos, -A_neg
end

A = matrix(ZZ, [[1, 2, 3], [4, 5, 6]])
r, n = size(A)


R, lambda = polynomial_ring(QQ, ["mu", "nu"])

function spower(lambda, a)
    r, n = size(a)
    @assert n == 1

    result = 1
    for i in 1:r
        result *= lambda[i]^(a[i, 1])
    end
    return result
end

function spower_matrix(lambda, A)
    r, n = size(A)

    result = [spower(lambda, A[:, i]) for i in 1:n]
end

function spower_neg(lambda, a)
    r, n = size(a)
    @assert n == 1
    a_pos, a_neg = split_matrix_pos_and_neg(a)
    res_pos = 1
    res_neg = 1
    for i in 1:r
        res_pos *= lambda[i]^(a_pos[i, 1])
        res_neg *= lambda[i]^(a_neg[i, 1])
    end
    return res_pos//res_neg
end

function spower_neg_matrix(lambda, A)
    r, n = size(A)

    result = [spower_neg(lambda, A[:, i]) for i in 1:n]
end

A = matrix(ZZ, [[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])

spower_neg_matrix(lambda, A)