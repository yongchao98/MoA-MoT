import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 0 or c <= 0:
        print("Error: n and c must be positive integers.")
        return

    # M1 is the number of mistakes when the true expert is wrong.
    # The true expert makes strictly fewer than c mistakes, so at most c-1.
    m1_bound = c - 1
    
    # M2 is the number of mistakes when the true expert is right.
    # For each such mistake, at least 2 bad experts must be wrong.
    # Total mistakes available to the (n-1) bad experts is (n-1)*c.
    # So, 2 * M2 <= (n-1)*c
    m2_bound = (n - 1) * c / 2
    
    # The total mistake bound is the sum of the bounds for M1 and M2.
    total_bound = m1_bound + m2_bound

    # To show the numbers in the final equation, we construct a formatted string.
    # We display both the formula with variables and the calculation with numbers.
    print(f"Given n = {n} experts and mistake threshold c = {c}")
    print("\nThe formula for the upper bound is: (c - 1) + ((n - 1) * c) / 2")
    
    print("\nCalculation:")
    # Printing each part of the equation with its calculated value
    term1_calc = f"({c} - 1)"
    term2_calc = f"(({n} - 1) * {c}) / 2"
    
    sub_term1_val = n - 1
    sub_term2_val = sub_term1_val * c
    
    print(f"{term1_calc} + {term2_calc}")
    print(f"= {m1_bound} + ({sub_term1_val} * {c}) / 2")
    print(f"= {m1_bound} + {sub_term2_val} / 2")
    print(f"= {m1_bound} + {m2_bound}")
    print(f"= {total_bound}")

    print(f"\nThe upper bound on the number of mistakes is: {total_bound}")


# --- Example Usage ---
# You can change these values to see the result for different scenarios.
num_experts = 11
mistake_threshold = 5

calculate_mistake_bound(num_experts, mistake_threshold)