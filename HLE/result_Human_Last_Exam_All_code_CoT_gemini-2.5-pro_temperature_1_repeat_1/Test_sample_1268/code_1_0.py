import math

def analyze_imaginary_quadratic_field(N):
    """
    Calculates and demonstrates the relationship between the Minkowski bound (M_k)
    and the covolume (V) for the imaginary quadratic field Q(sqrt(-N)).

    Args:
        N (int): A positive squarefree natural number.
    """
    print(f"--- Analysis for N = {N} ---")

    # Calculate the discriminant of the field K = Q(sqrt(-N))
    # If N = 1, 2 (mod 4), the discriminant is -4N.
    # If N = 3 (mod 4), the discriminant is -N.
    if N % 4 == 1 or N % 4 == 2:
        discriminant = -4 * N
    else:  # N % 4 == 3
        discriminant = -N

    abs_discriminant = abs(discriminant)
    
    # Calculate the covolume (V)
    # V = (1/2) * sqrt(|discriminant|)
    covolume = 0.5 * math.sqrt(abs_discriminant)

    # Calculate the Minkowski bound (M_k)
    # M_k = (2/pi) * sqrt(|discriminant|)
    minkowski_bound = (2 / math.pi) * math.sqrt(abs_discriminant)

    # The constant factor in the relationship M_k = C * V is C = 4/pi
    constant_factor = 4 / math.pi

    print(f"The imaginary quadratic field is K = Q(sqrt({-N})).")
    print(f"Field Discriminant (Δ_K): {discriminant}")
    print(f"Covolume (V): {covolume:.4f}")
    print(f"Minkowski Bound (k_k,∞): {minkowski_bound:.4f}")
    print("\nDerivation of the relationship:")
    print("The upper bound 'k_k,∞' is related to the covolume 'V' by the equation:")
    print("k_k,∞ = (4/π) * V")
    
    print("\nVerifying the equation with the calculated values:")
    # Print the equation with all the numbers
    print(f"{minkowski_bound:.4f} = {constant_factor:.4f} * {covolume:.4f}")
    
    # Calculate the right-hand side to show they are equal
    rhs = constant_factor * covolume
    print(f"{minkowski_bound:.4f} = {rhs:.4f}")
    
    # Final check
    if math.isclose(minkowski_bound, rhs):
        print("The relationship holds true.")
    else:
        print("There is a discrepancy in the calculation.")
    print("-" * 30 + "\n")


# --- Main execution ---
# We can test this for a few squarefree natural numbers.
# Example 1: N=3 (where discriminant is -N)
analyze_imaginary_quadratic_field(3)

# Example 2: N=5 (where discriminant is -4N)
analyze_imaginary_quadratic_field(5)

# Example 3: N=15 (composite but squarefree, 15 = 3*5)
analyze_imaginary_quadratic_field(15)
