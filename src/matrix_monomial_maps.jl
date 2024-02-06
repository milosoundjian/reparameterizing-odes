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

R, lambda = polynomial_ring(QQ, r, "lambda")

r, n = 3, 4

A = matrix(ZZ, [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])

apowermatrix(lambda, A)

function apower(lambda, a)
    r, n = size(a)

    result = one(lambda[1])
    
    for i in 1:r
        result *= lambda[i]^a[i]
    end

    return result
end

function apowermatrix(lambda, A)
    r, n = size(A)

    result = [one(lambda[1]) for i in 1:n]

    for i in 1:n
        result[i] = lmaopower(lambda, A[:, i])
    end

    return result
end

apower(lambda, a)

apowermatrix(lambda, A)

n

# a is a column vector of integers
a = matrix(ZZ, [[1], [2], [3]])
# lambda is a row vector with entries in K*
lambda = matrix(ZZ, [[2, 3, 4]])
spower(lambda, a)
