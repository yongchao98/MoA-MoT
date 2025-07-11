import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes made by a
    majority voting algorithm with expert removal.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(c, int) or c <= 0:
        print("Error: c must be a positive integer.")
        return

    # M1 is the number of mistakes where the true expert is also wrong.
    # The true expert makes strictly fewer than c mistakes.
    m1_bound = c - 1

    # M2 is the number of mistakes where the true expert is right.
    # In this case, at least two "false" experts must be wrong.
    # The total mistake budget for n-1 false experts is (n-1)*c.
    # So, 2 * M2 <= (n-1)*c
    m2_bound = (n - 1) * c // 2  # Using floor division for integer result

    total_bound = m1_bound + m2_bound

    print("An upper bound on the number of algorithm mistakes is the sum of two cases:")
    print("1. M1: Algorithm and true expert are both wrong.")
    print("2. M2: Algorithm is wrong, but the true expert is right.")
    print("-" * 20)
    
    # Print the breakdown of the calculation as requested
    print(f"Given n = {n} experts and mistake threshold c = {c}:")
    print("\nCalculating the bound for M1:")
    print(f"M1 <= c - 1")
    print(f"M1 <= {c} - 1 = {m1_bound}")

    print("\nCalculating the bound for M2:")
    print(f"M2 <= floor((n - 1) * c / 2)")
    print(f"M2 <= floor(({n} - 1) * {c} / 2) = {m2_bound}")

    print("\nCalculating the total upper bound:")
    print(f"Total Mistakes <= M1 + M2")
    print(f"Total Mistakes <= {m1_bound} + {m2_bound} = {total_bound}")


# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
num_experts = 10
mistake_threshold = 5
calculate_mistake_bound(num_experts, mistake_threshold)