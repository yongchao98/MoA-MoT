import math

def calculate_mistake_bound(n, c):
    """
    Calculates and explains the upper bound on the number of mistakes for the given
    experts problem variant.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be integers greater than or equal to 1.")
        return

    print(f"Calculating the upper bound on algorithm mistakes for n={n} experts and a mistake limit of c={c}.")
    print("-" * 50)

    # Bound for M1 (mistakes when true expert is wrong)
    m1_bound = c - 1
    print("Step 1: Bound mistakes when the true expert is also wrong (M1).")
    print(f"The true expert makes at most c-1 mistakes.")
    print(f"M1 <= c - 1")
    print(f"M1 <= {c} - 1 = {m1_bound}")
    print("-" * 50)

    # Bound for M2 (mistakes when true expert is right)
    non_true_experts = n - 1
    total_non_true_mistakes = non_true_experts * c
    m2_bound = math.floor(total_non_true_mistakes / 2) # or use integer division //

    print("Step 2: Bound mistakes when the true expert is right (M2).")
    print(f"The total number of mistakes the other {non_true_experts} experts can make is at most ({n} - 1) * {c} = {total_non_true_mistakes}.")
    print(f"For each M2 mistake, the number of wrong non-true experts is at least 2 more than the right ones.")
    print(f"This implies that 2 * M2 is less than or equal to the total mistakes made by non-true experts.")
    print(f"2 * M2 <= ({n} - 1) * {c}")
    print(f"M2 <= ({non_true_experts} * {c}) / 2")
    print(f"M2 <= {total_non_true_mistakes} / 2 = {total_non_true_mistakes / 2}")
    print(f"Since M2 is an integer, M2 <= floor({total_non_true_mistakes / 2}) = {m2_bound}")
    print("-" * 50)

    # Total Bound
    total_bound = m1_bound + m2_bound
    print("Step 3: Combine the bounds for the total number of mistakes (M).")
    print("M <= M1 + M2")
    print(f"M <= {m1_bound} + {m2_bound}")
    print(f"M <= {total_bound}")
    print("-" * 50)
    print(f"The calculated upper bound on the number of mistakes is: {total_bound}")


# --- User Input ---
# You can change these values to test different scenarios
n_experts = 11
mistake_limit = 4

calculate_mistake_bound(n_experts, mistake_limit)
