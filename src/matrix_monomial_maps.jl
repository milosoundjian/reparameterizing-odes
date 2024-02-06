using Nemo

function spower(lambda, a)
    # Initialize the result
    result = one(lambda[1])

    # Compute the product of lambda_i^a_i
    for i in axes(lambda, 2)
        result *= lambda[1, i]^a[i, 1]
    end

    return result
end

function spower_matrix(lambda, A)
    r, n = size(A)
    result = matrix(ZZ, [[one(lambda[1]) for i in 1:n]])

    for i in 1:n
        a = A[:, i]
        result[i] = spower(lambda, a)
    end

    return result
end

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

r, n = 2, 3
A = matrix(ZZ, [[1, 2, 3], [4, 5, 6]])

lambda = matrix(ZZ, [[2, 3]])

spower_matrix(lambda, A)

A_pos, A_neg = split_matrix_pos_and_neg(A)