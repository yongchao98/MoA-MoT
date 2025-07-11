import math

def get_upper_bound_relation(N: int):
    """
    Calculates and displays the relationship between the maximum norm and covolume
    for the real quadratic field K = Q(sqrt(N)).

    Args:
        N (int): A squarefree natural number.
    """
    # 1. Validate that N is a squarefree natural number.
    if not isinstance(N, int) or N <= 0:
        print(f"Error: Input '{N}' must be a natural number (positive integer).")
        return

    i = 2
    temp_n = N
    while i * i <= temp_n:
        if temp_n % (i * i) == 0:
            print(f"Error: Input N={N} is not squarefree as it is divisible by {i*i}.")
            return
        i += 1

    print(f"Analyzing the field K = Q(sqrt({N}))...")

    # 2. Calculate the discriminant D_K of the number field.
    if N % 4 == 1:
        discriminant = N
    else:  # This covers N % 4 == 2 and N % 4 == 3
        discriminant = 4 * N
    print(f"The discriminant of the field is D_K = {discriminant}.")

    # 3. Calculate the covolume V.
    covolume = math.sqrt(discriminant)
    print(f"The covolume is V = sqrt({discriminant}) ≈ {covolume:.4f}.")

    # 4. The Minkowski bound gives the relation k_k,inf <= (1/2) * V.
    constant = 0.5
    bound_value = constant * covolume

    print(f"\nThe upper bound for the maximum norm (k_k,∞) is {constant} times the covolume V.")
    print(f"For K = Q(sqrt({N})), this value is {constant} * {covolume:.4f} ≈ {bound_value:.4f}.")

    # 5. Output the final relationship, showing each number in the equation.
    print("\nThe general relationship is expressed by the inequality:")
    print(f"k_k,∞ <= {constant} * V")


# --- Main execution ---
# We will use N=15 as an example squarefree natural number.
# 15 = 3 * 5, so it is squarefree.
# 15 % 4 = 3, so the discriminant D_K will be 4 * 15 = 60.
example_N = 15
get_upper_bound_relation(example_N)
