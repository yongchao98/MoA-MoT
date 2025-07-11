import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by the algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert makes before being removed.
    """
    if n <= 1 or c <= 0:
        print("n must be > 1 and c must be > 0")
        return

    # Bound for mistakes where the true expert was right
    # M_good <= floor(c * (n - 1) / 2)
    m_good_bound = math.floor(c * (n - 1) / 2)

    # Bound for mistakes where the true expert was wrong
    # M_bad <= c - 1
    m_bad_bound = c - 1

    # Total bound is the sum of the two bounds
    total_bound = m_good_bound + m_bad_bound

    print("Derivation of the upper bound for the algorithm's mistakes:")
    print(f"For n = {n} experts and a mistake cutoff c = {c}:")
    print("\n1. Bound on mistakes when the true expert is right (M_good):")
    print(f"   M_good <= floor(c * (n - 1) / 2)")
    print(f"   M_good <= floor({c} * ({n} - 1) / 2) = floor({c * (n-1)} / 2) = {m_good_bound}")

    print("\n2. Bound on mistakes when the true expert is wrong (M_bad):")
    print(f"   M_bad <= c - 1")
    print(f"   M_bad <= {c} - 1 = {m_bad_bound}")

    print("\n3. Total Upper Bound (M):")
    print(f"   M = M_good + M_bad")
    print(f"   M <= {m_good_bound} + {m_bad_bound} = {total_bound}")
    print(f"\nThe final upper bound is {total_bound}.")


# Example values for n and c
n_experts = 11
mistake_cutoff = 10

calculate_mistake_bound(n_experts, mistake_cutoff)
