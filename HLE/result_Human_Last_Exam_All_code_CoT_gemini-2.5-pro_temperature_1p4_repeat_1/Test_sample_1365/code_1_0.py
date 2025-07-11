import math

def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes for the given expert problem variant.

    Args:
        n (int): The number of experts.
        c (int): The mistake threshold for removing an expert.
    """

    print(f"Calculating the upper bound for n = {n} experts and a mistake threshold c = {c}.\n")

    if n <= 0 or c <= 0:
        print("Number of experts (n) and mistake threshold (c) must be positive integers.")
        return

    # M2 is the number of mistakes where the true expert is also wrong.
    # The true expert makes strictly fewer than c mistakes.
    m2_bound = c - 1
    print(f"Step 1: Bound on M2 (mistakes when the true expert is wrong)")
    print(f"M2 <= c - 1")
    print(f"M2 <= {c} - 1 = {m2_bound}\n")


    # M1 is the number of mistakes where the true expert is correct.
    # This can only happen if other experts outvote the true expert.
    # For a strict majority vote, this requires at least 2 other experts to be wrong.
    # This means M1 can only be non-zero if n > 2.
    if n <= 2:
        m1_bound = 0
        print(f"Step 2: Bound on M1 (mistakes when the true expert is correct)")
        print(f"With n = {n}, the true expert cannot be outvoted by a strict majority.")
        print(f"Therefore, M1 = {m1_bound}\n")
    else:
        # Total mistakes available to the n-1 "other" experts is (n-1)*c.
        # Each M1-type mistake requires at least 2 mistakes from this pool.
        m1_bound = math.floor((n - 1) * c / 2)
        print(f"Step 2: Bound on M1 (mistakes when the true expert is correct)")
        print(f"Each M1 mistake requires at least 2 'other' experts to be wrong.")
        print(f"The mistake pool for the {n-1} other experts is ({n} - 1) * {c} = {(n-1)*c}.")
        print(f"M1 <= floor( (n - 1) * c / 2 )")
        print(f"M1 <= floor( ({n} - 1) * {c} / 2 ) = floor({(n - 1) * c / 2}) = {m1_bound}\n")

    # The total mistake bound is M1 + M2
    total_bound = m1_bound + m2_bound
    print(f"Step 3: Total mistake bound M = M1 + M2")
    print(f"M <= {m1_bound} + {m2_bound} = {total_bound}\n")
    print(f"The final upper bound on the number of mistakes is: {total_bound}")

# Example usage:
# You can change these values to test with different n and c
n_experts = 10
c_threshold = 5

calculate_mistake_bound(n_experts, c_threshold)
