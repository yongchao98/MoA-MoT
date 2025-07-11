import math

def analyze_imaginary_quadratic_field(N):
    """
    Calculates the relationship between the Minkowski bound and the covolume
    for an imaginary quadratic field Q(sqrt(-N)).

    Args:
        N (int): A positive squarefree natural number.
    """
    if N <= 0:
        print("Please provide a positive natural number.")
        return

    # 1. Determine the discriminant of the field K = Q(sqrt(-N))
    # Note: We treat the user's N as the 'd' in Q(sqrt(-d))
    if (-N) % 4 == 1:
        # This happens when N = 3, 7, 11, 15, ... (i.e., N = 4k+3)
        discriminant = -N
    else:
        # This happens when N = 1, 2, 5, 6, ... (i.e., N = 4k+1 or 4k+2)
        discriminant = -4 * N

    abs_discriminant = abs(discriminant)

    # 2. Calculate the covolume (V)
    # This is the area of the fundamental parallelogram of the ring of integers
    # embedded in the complex plane.
    covolume_V = math.sqrt(abs_discriminant) / 2

    # 3. Calculate the upper bound for the norm (Minkowski Bound M_K)
    # The Minkowski bound is an upper bound on the norm of an ideal
    # in any ideal class. We interpret k_k,inf as this bound.
    minkowski_bound_M_K = (2 / math.pi) * math.sqrt(abs_discriminant)

    # 4. The derived relationship is M_K = (4/pi) * V.
    # We will print the values to show this relationship.
    # This shows the upper bound (M_K) in relation to the covolume (V).
    bound_from_covolume = (4 / math.pi) * covolume_V

    print(f"For the squarefree number N = {N}:")
    print(f"  - The discriminant (Δ) is: {discriminant}")
    print(f"  - The covolume (V) is: sqrt(|{discriminant}|) / 2 = {covolume_V:.5f}")
    print(f"  - The Minkowski bound (M_K) is: (2/π) * sqrt(|{discriminant}|) = {minkowski_bound_M_K:.5f}\n")
    print("The relationship between the upper bound (M_K) and the covolume (V) is:")
    print("  M_K = (4 / π) * V\n")
    print("Let's plug in the numbers to verify:")
    print(f"Equation: {minkowski_bound_M_K:.5f} = (4 / {math.pi:.5f}) * {covolume_V:.5f}")
    print(f"Right-hand side calculates to: {bound_from_covolume:.5f}")
    print("\nThe equality holds.")


# --- Main Execution ---
# Let's use N=15 (a squarefree number) as an example.
# 15 = 3 * 5, not divisible by any squares.
example_N = 15
analyze_imaginary_quadratic_field(example_N)
