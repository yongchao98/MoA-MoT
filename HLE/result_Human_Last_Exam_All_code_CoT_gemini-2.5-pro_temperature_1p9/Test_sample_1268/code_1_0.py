import math

def calculate_minkowski_bound(N):
    """
    Calculates the Minkowski bound for the real quadratic field Q(sqrt(N)).
    
    Args:
    N (int): A positive squarefree integer.
    
    Returns:
    None. Prints the step-by-step calculation.
    """
    print(f"Starting with the squarefree natural number N = {N}")
    print("-" * 30)

    # Step 1: Calculate the discriminant of the field Q(sqrt(N))
    # The discriminant is N if N = 1 (mod 4)
    # The discriminant is 4N if N = 2 or 3 (mod 4)
    if N % 4 == 1:
        discriminant = N
    else:
        discriminant = 4 * N
    print(f"Step 1: The discriminant (Δ) for Q(sqrt({N})) is {discriminant}.")

    # Step 2: Calculate the covolume V
    # For a real quadratic field, the covolume V is the square root of the discriminant.
    covolume_V = math.sqrt(discriminant)
    print(f"Step 2: The covolume (V) is √|Δ| = √{discriminant} ≈ {covolume_V:.4f}.")

    # Step 3: State the relationship and calculate the bound
    # The upper bound for the norm of an ideal in any class is (1/2) * V
    # We interpret k_{k,∞} to be this bound.
    minkowski_constant = 0.5
    upper_bound = minkowski_constant * covolume_V
    print(f"Step 3: The Minkowski bound is calculated as {minkowski_constant} * V.")
    
    print("\n" + "=" * 30)
    print("Final Result")
    print("=" * 30)
    print("The upper bound for the norm, k_k,∞, is related to the covolume, V, by the inequality:")
    print(f"k_k,∞ <= C * V")
    print(f"where C is the constant from the Minkowski bound.")
    print("\nIn this case, for a real quadratic field:")
    print(f"k_k,∞ <= {minkowski_constant} * V")
    
    print("\nFor our specific example with N=15:")
    print(f"The equation becomes k_k,∞ <= {minkowski_constant} * {covolume_V:.4f}")
    print(f"Therefore, the upper bound is approximately {upper_bound:.4f}")


# Use a squarefree number as an example, e.g., N=15 (3*5)
squarefree_N = 15
calculate_minkowski_bound(squarefree_N)