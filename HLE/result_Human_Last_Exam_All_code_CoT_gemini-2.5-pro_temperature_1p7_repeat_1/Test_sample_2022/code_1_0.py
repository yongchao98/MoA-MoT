def mat_mul(A, B, p):
    """Multiplies two 3x3 matrices A and B modulo p."""
    C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % p
    return C

def mat_pow(A, n, p):
    """Computes A^n modulo p for a 3x3 matrix A using exponentiation by squaring."""
    res = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        res[i][i] = 1  # Identity matrix
    base = A
    while n > 0:
        if n % 2 == 1:
            res = mat_mul(res, base, p)
        base = mat_mul(base, base, p)
        n //= 2
    return res

def solve():
    # Case p = 80039
    p1 = 80039
    # The value is S(3) = 7/8 mod p1
    # We calculate 7 * (8^-1 mod p1)
    inv8_p1 = pow(8, p1 - 2, p1)
    ans1 = (7 * inv8_p1) % p1

    # Case p = 80077
    p2 = 80077
    # The value is S((p2-1)^2) mod p2
    # Recurrence: S_k = S_{k-1} + c * S_{k-3}, with c = -1/8
    # Matrix form: v_k = M * v_{k-1}, where v_k = [S_k, S_{k-1}, S_{k-2}]^T
    # M = [[1, 0, c], [1, 0, 0], [0, 1, 0]]
    # We need to calculate v_k = M^(k-2) * v_2
    # where k = (p2-1)^2 and v_2 = [S_2, S_1, S_0]^T = [1, 1, 1]^T
    c = (-pow(8, p2 - 2, p2)) % p2
    M = [[1, 0, c], [1, 0, 0], [0, 1, 0]]
    k = (p2 - 1)**2
    exp = k - 2
    
    if exp < 0:
        # This case won't be hit with the current problem values, but good practice
        # M_inv = ...
        # M_E = mat_pow(M_inv, -exp, p2)
        # Not needed here.
        pass
    else:
        M_E = mat_pow(M, exp, p2)
    
    # S_k is the first component of M_E * [1, 1, 1]^T
    # S_k = M_E[0][0]*S2 + M_E[0][1]*S1 + M_E[0][2]*S0
    # S_k = M_E[0][0] + M_E[0][1] + M_E[0][2]
    ans2 = (M_E[0][0] + M_E[0][1] + M_E[0][2]) % p2

    print(f"{ans1},{ans2}")

solve()