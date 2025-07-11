import math

def analyze_quadratic_field_bounds(N):
    """
    Analyzes and demonstrates the relationship between the Minkowski bound (M_K) 
    and the covolume (V) for quadratic fields derived from a squarefree number N.
    
    The term "maximum norm (k_k,∞)" is interpreted as the Minkowski Bound.
    """
    print(f"Analyzing for the squarefree number N = {N}...\n")

    # --- Case 1: Real Quadratic Field K = Q(sqrt(N)) ---
    print("--- Case 1: Real Quadratic Field K = Q(sqrt(N)) ---")
    d_real = N
    # Calculate discriminant Δ_K for the real field
    if d_real % 4 == 1:
        delta_real = d_real
    else:
        delta_real = 4 * d_real
    
    n = 2      # Degree of the field
    r2_real = 0  # Number of pairs of complex embeddings
    
    # Calculate Covolume V = sqrt(|Δ_K|) / 2^r2
    V_real = math.sqrt(delta_real)
    
    # Calculate Minkowski Bound M_K = (n!/n^n) * (4/π)^r2 * sqrt(|Δ_K|)
    M_K_real = (math.factorial(n) / (n**n)) * math.sqrt(delta_real)

    # The constant C in the relation M_K = C * V
    C_real = 0.5

    print(f"For K = Q(sqrt({N})), the field parameters are:")
    print(f"  - Discriminant Δ_K = {delta_real}")
    print(f"  - Covolume V = sqrt({delta_real}) ≈ {V_real:.5f}")
    print(f"  - Minkowski Bound M_K = (1/2)*sqrt({delta_real}) ≈ {M_K_real:.5f}")
    print("  The relationship is M_K = C * V, where C = 1/2.")
    print(f"  Equation: {M_K_real:.5f} = {C_real} * {V_real:.5f}")
    print("-" * 50)

    # --- Case 2: Imaginary Quadratic Field K = Q(sqrt(-N)) ---
    print("\n--- Case 2: Imaginary Quadratic Field K = Q(sqrt(-N)) ---")
    d_imag = -N
    # Calculate discriminant Δ_K for the imaginary field
    if d_imag % 4 == 1:
        delta_imag = d_imag
    else:
        delta_imag = 4 * d_imag
    
    delta_imag_abs = abs(delta_imag)
    r2_imag = 1  # Number of pairs of complex embeddings

    # Calculate Covolume V = sqrt(|Δ_K|) / 2^r2
    V_imag = math.sqrt(delta_imag_abs) / 2.0
    
    # Calculate Minkowski Bound M_K = (n!/n^n) * (4/π)^r2 * sqrt(|Δ_K|)
    M_K_imag = (math.factorial(n) / (n**n)) * (4 / math.pi) * math.sqrt(delta_imag_abs)

    # The constant C in the relation M_K = C * V
    C_imag = 4 / math.pi
    
    print(f"For K = Q(sqrt({d_imag})), the field parameters are:")
    print(f"  - Discriminant Δ_K = {delta_imag}")
    print(f"  - Covolume V = sqrt({delta_imag_abs})/2 ≈ {V_imag:.5f}")
    print(f"  - Minkowski Bound M_K = (2/π)*sqrt({delta_imag_abs}) ≈ {M_K_imag:.5f}")
    print(f"  The relationship is M_K = C * V, where C = 4/π ≈ {C_imag:.5f}")
    print(f"  Equation: {M_K_imag:.5f} = {C_imag:.5f} * {V_imag:.5f}")
    print("-" * 50)

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print("The general upper bound is determined by the larger of the two constants (1/2 and 4/π).")
    print(f"Since 4/π ≈ {C_imag:.5f}, which is greater than 0.5, this value sets the bound.")
    print("The upper bound for the Minkowski norm (M_K) in relation to the covolume (V) is:")
    print("  M_K ≤ (4/π) * V")
    print(f"\nThe constant factor defining this upper bound is 4/π ≈ {4/math.pi:.5f}.")


def is_squarefree(n):
    """Helper function to check if a number is squarefree."""
    if n <= 0: return False
    i = 2
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True

# --- Main execution ---
# We use a squarefree number like 105 (3 * 5 * 7) as an example.
example_N = 105
if is_squarefree(example_N):
    analyze_quadratic_field_bounds(example_N)
else:
    print(f"Error: The number {example_N} is not squarefree. Please use a squarefree number.")
