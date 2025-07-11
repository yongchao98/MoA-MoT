def solve():
    """
    Solves the coloring problem by using a linear recurrence relation and modular arithmetic.
    """
    # Prime modulus
    p = 23627
    print(f"The modulus p is {p}, which is a prime number.")

    # Number of colors
    num_colors = 510
    num_bad_colors = 203

    # M is the number of ways to color a single 2x1 column
    M = num_colors**2
    # A is the number of "bad" monochromatic columns
    A = num_bad_colors

    M_mod_p = M % p
    A_mod_p = A % p
    print(f"Total column colorings M = {num_colors}^2 = {M}. M mod {p} = {M_mod_p}.")
    print(f"Bad mono-column types A = {num_bad_colors}. A mod {p} = {A_mod_p}.")

    # The recurrence is S(n) = (M-1)S(n-1) + (M-1)S(n-2) + (M-A)S(n-3)
    # Modulo p, M-A = 0.
    # So, S(n) = (M-1)S(n-1) + (M-1)S(n-2) mod p for n >= 3.
    c = (M_mod_p - 1 + p) % p
    print(f"The recurrence becomes S(n) = {c}(S(n-1) + S(n-2)) mod {p} for n>=3.")

    # Calculate initial terms of the sequence
    S0 = 1
    S1 = M_mod_p
    S2 = (M_mod_p**2) % p
    # S3 is the first term governed by the simplified recurrence
    S3 = (c * (S2 + S1)) % p
    
    print(f"Initial values mod {p}:")
    print(f"S(0) = {S0}")
    print(f"S(1) = {S1}")
    print(f"S(2) = {S2}")
    print(f"S(3) = {S3}")

    # The index n = 23626 * (23628^100 - 23628^50) is a multiple of the
    # sequence period. We need to compute C^{n-3}, which simplifies to C^{-3}.
    # The recurrence matrix C = [[c, c], [1, 0]].
    # We compute C inverse.
    c_inv = pow(c, -1, p)
    
    # C_inv = [[0, 1], [c_inv, -1]]
    C_inv = [[0, 1], [c_inv, (p - 1)]]

    # Calculate C^{-2} = C^{-1} * C^{-1}
    C_inv2_11 = (C_inv[0][0] * C_inv[0][0] + C_inv[0][1] * C_inv[1][0]) % p
    C_inv2_12 = (C_inv[0][0] * C_inv[0][1] + C_inv[0][1] * C_inv[1][1]) % p
    C_inv2_21 = (C_inv[1][0] * C_inv[0][0] + C_inv[1][1] * C_inv[1][0]) % p
    C_inv2_22 = (C_inv[1][0] * C_inv[0][1] + C_inv[1][1] * C_inv[1][1]) % p
    C_inv2 = [[C_inv2_11, C_inv2_12], [C_inv2_21, C_inv2_22]]

    # Calculate C^{-3} = C^{-1} * C^{-2}
    C_inv3_11 = (C_inv[0][0] * C_inv2[0][0] + C_inv[0][1] * C_inv2[1][0]) % p
    C_inv3_12 = (C_inv[0][0] * C_inv2[0][1] + C_inv[0][1] * C_inv2[1][1]) % p
    
    print(f"The state vector evolves according to V_n = C^(n-3) * V_3 for n>=3.")
    print(f"The index n is a multiple of the period of C, so we compute C^(-3).")
    print(f"First row of C^(-3) mod {p} is: ({C_inv3_11}, {C_inv3_12}).")

    # S(n) is the first component of V_n = C^{n-3} * V_3, where V_3 = [S3, S2]^T
    Sn = (C_inv3_11 * S3 + C_inv3_12 * S2) % p

    n_str = "23626*(23628**100-23628**50)"
    p_val = 23627

    # Final result
    print("\nFinal calculation:")
    print(f"S(n) mod p = (C^(-3) * [S(3), S(2)]^T)_1")
    print(f"S(n) mod {p_val} = ({C_inv3_11} * {S3} + {C_inv3_12} * {S2}) mod {p_val}")
    print(f"S(n) mod {p_val} = ({(C_inv3_11 * S3) % p_val} + {(C_inv3_12 * S2) % p_val}) mod {p_val}")
    print(f"S(n) mod {p_val} = {Sn}")
    print("\nFinal equation:")
    print(f"S({n_str}) mod {p_val} = {Sn}")
    return Sn

final_answer = solve()
<<<11613>>>