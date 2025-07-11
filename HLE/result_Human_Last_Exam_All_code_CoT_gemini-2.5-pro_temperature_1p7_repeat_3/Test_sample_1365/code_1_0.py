import math

def calculate_mistake_bound(n, c):
    """
    Calculates and prints an upper bound on the number of mistakes for the experts problem variant.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("n and c must be integers greater than or equal to 1.")
        return

    # Bound for Type 1 mistakes (M1), where the true expert is wrong.
    # This is limited by the number of mistakes the true expert can make.
    m1_bound = c - 1

    # Bound for Type 2 mistakes (M2), where the true expert is right.
    # This is limited by the total mistake budget of the n-1 false experts.
    m2_bound = (c * (n - 1)) / 2

    # Total bound is the sum of the two bounds.
    total_bound = m1_bound + m2_bound

    print("Derivation of the upper bound on algorithm mistakes:")
    print("-" * 50)
    print(f"Given n = {n} experts and a mistake limit c = {c}.")
    print("\nStep 1: Bound on mistakes when the true expert is also wrong (M1).")
    print(f"The true expert makes at most c-1 mistakes.")
    print(f"M1 <= {c} - 1 = {m1_bound}")

    print("\nStep 2: Bound on mistakes when the true expert is right (M2).")
    print(f"The total mistake budget for the {n-1} false experts is {c} * ({n} - 1) = {c * (n - 1)}.")
    print("For each such mistake, the number of wrong false experts exceeds right false experts by at least 2.")
    print(f"This leads to M2 <= (c * (n - 1)) / 2.")
    print(f"M2 <= ({c} * ({n} - 1)) / 2 = {m2_bound}")

    print("\nStep 3: Total Upper Bound (M).")
    print(f"M <= M1 + M2")
    # Using floor for the final bound since the number of mistakes is an integer.
    # The calculated value is an upper bound, so floor is appropriate.
    final_bound_int = math.floor(total_bound)
    print(f"M <= ({c} - 1) + ({c} * ({n} - 1)) / 2 = {m1_bound} + {m2_bound} = {total_bound}")
    print(f"Since the number of mistakes must be an integer, the upper bound is {final_bound_int}.")


# Example usage with n=11 experts and c=4 mistake threshold.
calculate_mistake_bound(n=11, c=4)