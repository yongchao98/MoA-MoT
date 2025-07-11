import math

# --- Parameters of the problem ---
# n: the total number of experts
n = 11

# c: the number of mistakes after which an expert is removed
c = 5
# The true expert makes strictly fewer than c mistakes.
# ------------------------------------

def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes
    made by the majority voting algorithm.
    """
    if n <= 0 or c <= 0:
        print("n and c must be positive integers.")
        return

    # Bound for M_0 (mistakes when the true expert is wrong)
    # M_0 <= c - 1
    m0_bound = c - 1

    # Bound for M_1 (mistakes when the true expert is right)
    # M_1 <= floor(c * (n - 1) / 2)
    m1_bound_float = c * (n - 1) / 2
    m1_bound_int = math.floor(m1_bound_float)

    # Total bound M = M_0 + M_1
    total_bound = m0_bound + m1_bound_int

    print("An upper bound on the number of mistakes (M) is given by the formula:")
    print("M <= (c - 1) + floor(c * (n - 1) / 2)")
    print("\nFor the given values:")
    print(f"n = {n}")
    print(f"c = {c}")
    print("\nThe calculation is as follows:")
    print(f"M <= ({c} - 1) + floor({c} * ({n} - 1) / 2)")
    print(f"M <= {m0_bound} + floor({c * (n - 1)} / 2)")
    print(f"M <= {m0_bound} + floor({m1_bound_float})")
    print(f"M <= {m0_bound} + {m1_bound_int}")
    print(f"M <= {total_bound}")
    print("\nTherefore, the algorithm will make at most", total_bound, "mistakes.")

calculate_mistake_bound(n, c)