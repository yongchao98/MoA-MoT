def mat_mul(A, B):
    """Multiplies two matrices A and B."""
    C = [[0, 0], [0, 0]]
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
    return C

def mat_pow(A, n):
    """Computes the n-th power of a 2x2 matrix A using exponentiation by squaring."""
    res = [[1, 0], [0, 1]]  # Identity matrix
    base = A
    
    while n > 0:
        if n % 2 == 1:
            res = mat_mul(res, base)
        base = mat_mul(base, base)
        n //= 2
    return res

def mat_vec_mul(A, v):
    """Multiplies a 2x2 matrix A by a 2x1 vector v."""
    res = [0, 0]
    res[0] = A[0][0] * v[0] + A[0][1] * v[1]
    res[1] = A[1][0] * v[0] + A[1][1] * v[1]
    return res

# The matrix M from the recurrence relation S_N = 4*S_{N-1} + 2*T_{N-1} and T_N = 3*S_{N-1} + 2*T_{N-1}
M = [[4, 2], [3, 2]]

# The initial vector v_0 = [S_0, T_0]
v0 = [4, 3]

# We need to find S_19, which corresponds to N=19
power = 19

# Compute M^19
M_pow_19 = mat_pow(M, power)

# Compute the final vector v_19 = M^19 * v_0
v19 = mat_vec_mul(M_pow_19, v0)

# The result is the first component of the vector, S_19
result = v19[0]

# Print the final numerical answer
print(f"The recurrence relation for S_N (the sum of squares for P_N(x)) and T_N (the coefficient of x^1 in P_N(x)P_N(x^-1)) is:")
print(f"S_N = 4 * S_(N-1) + 2 * T_(N-1)")
print(f"T_N = 3 * S_(N-1) + 2 * T_(N-1)")
print(f"With initial values S_0 = {v0[0]} and T_0 = {v0[1]}.")
print(f"We need to compute S_19, which is the first component of M^19 * v_0, where M is the matrix of the recurrence.")
print(f"The value of M^19 is: {M_pow_19[0]}")
print(f"                      {M_pow_19[1]}")
print(f"The final result for the sum of squares, S_19, is:")
print(result)