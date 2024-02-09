using Nemo

# Verhulst Model of Logistic Growth
# dn/dt = r * n * (1 - n/k)

# Let us get the Laurent Polynomial
R, (r, k, t, n) = laurent_polynomial_ring(QQ, ["r", "k", "t", "n"])

F = [r * t, r * k^(-1) * t * n]

S, (r1, k1, t1, n1) = laurent_polynomial_ring(QQ, ["r1", "k1", "t1", "n1"])

v = [r1 * 0 + 1, k1, t1, n1]

result = [evaluate(F[i], v) for i in 1:2]