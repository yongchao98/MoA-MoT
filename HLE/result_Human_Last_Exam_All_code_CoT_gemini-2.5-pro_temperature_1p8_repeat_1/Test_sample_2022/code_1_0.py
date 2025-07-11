def matrix_mult(A, B, p):
    """Multiplies two 3x3 matrices modulo p."""
    C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % p
    return C

def matrix_pow(A, n, p):
    """Computes matrix A to the power of n modulo p."""
    if n == 0:
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Identity matrix
    if n < 0:
        # For negative exponent, compute power of inverse matrix
        det = 0
        for i in range(3):
            det = (det + (A[0][i] * (A[1][(i + 1) % 3] * A[2][(i + 2) % 3] - A[1][(i + 2) % 3] * A[2][(i + 1) % 3]))) % p
        
        inv_det = pow(det, -1, p)
        
        A_inv = [[0,0,0],[0,0,0],[0,0,0]]
        for i in range(3):
          for j in range(3):
            A_inv[j][i] = ( (A[(i+1)%3][(j+1)%3] * A[(i+2)%3][(j+2)%3]) - (A[(i+1)%3][(j+2)%3] * A[(i+2)%3][(j+1)%3]) ) * inv_det % p
        
        return matrix_pow(A_inv, -n, p)

    # Exponentiation by squaring
    res = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    base = A
    while n > 0:
        if n % 2 == 1:
            res = matrix_mult(res, base, p)
        base = matrix_mult(base, base, p)
        n //= 2
    return res

def solve():
    primes = [80039, 80077]
    results = []
    for p in primes:
        # The recurrence relation is F(n) = F(n-1) - 1/8 * F(n-3) mod p
        # Which is F(n+3) = F(n+2) - 1/8 * F(n) mod p
        inv8 = pow(8, -1, p)
        neg_inv8 = (-inv8 + p) % p

        # Companion matrix for state vector [F(n), F(n+1), F(n+2)]^T
        C = [
            [0, 1, 0],
            [0, 0, 1],
            [neg_inv8, 0, 1]
        ]

        if p == 80039:
            # For p=80039, 5 is a quadratic residue, period divides p-1.
            # N = p^5+2p^4-19p^3-3p^2+16p+6
            # N mod (p-1) = 1+2-19-3+16+6 = 3
            n = 3
        else: # p == 80077
            # For p=80077, 5 is a quadratic non-residue, period divides p^2-1.
            # N mod (p^2-1) = p(1) + 2(1) - 19p(1) - 3(1) + 16p + 6
            # N mod (p^2-1) = p + 2 - 19p - 3 + 16p + 6 = -2p + 5
            n = -2 * p + 5
        
        # We need to calculate F(n). We use matrix exponentiation
        # v_n = C^n * v_0, where v_0 = [F(0), F(1), F(2)] = [1, 1, 1]
        # F(n) is the first component of v_n
        if n == 3: # Short-circuit for small n
            # F(0)=1, F(1)=1, F(2)=1
            # F(3) = F(2) - inv8 * F(0) = 1 - inv8
            f_n = (1 - inv8 + p) % p
            # Equation: F(3) = 1 - 8^(-1) mod 80039
            val_inv8 = pow(8, -1, 80039)
            val_final = (1 - val_inv8 + 80039) % 80039
            print(f"For p=80039, the value is F(3) = (1 - {val_inv8}) mod 80039 = {val_final}")

        else:
            C_n = matrix_pow(C, n, p)
            # v_0 = [1, 1, 1]
            f_n = (C_n[0][0] * 1 + C_n[0][1] * 1 + C_n[0][2] * 1) % p
            
            # Print intermediate values
            exponent = -2 * 80077 + 5
            print(f"For p=80077, the value is F({exponent}) mod {80077} = {f_n}")

        results.append(f_n)

    print(f"\nFinal numerical answers separated by a comma:")
    print(f"{results[0]},{results[1]}")

solve()
<<<70035,50049>>>