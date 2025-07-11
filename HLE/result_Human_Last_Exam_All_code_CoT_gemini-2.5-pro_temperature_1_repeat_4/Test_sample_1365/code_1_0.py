import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for the described
    majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be integers greater than 0.")
        return

    # Part 1: Bound for M_wrong (mistakes when the true expert is wrong)
    # The true expert makes at most c-1 mistakes.
    m1_bound = c - 1

    # Part 2: Bound for M_correct (mistakes when the true expert is correct)
    # This bound is derived from the total mistake budget of the false experts.
    term_in_floor_num = (n - 1) * c
    term_in_floor = term_in_floor_num / 2
    m2_bound = math.floor(term_in_floor)

    # Total upper bound
    total_bound = m1_bound + m2_bound

    # --- Output ---
    print(f"For n = {n} experts and a mistake limit c = {c}:")
    print("The upper bound on the number of mistakes is given by the formula:")
    print("Bound = (c - 1) + floor((n - 1) * c / 2)")
    print("\nCalculation steps:")
    # Here we output each number in the final equation
    print(f"Bound = ({c} - 1) + floor(({n} - 1) * {c} / 2)")
    print(f"Bound = {m1_bound} + floor({term_in_floor_num} / 2)")
    print(f"Bound = {m1_bound} + floor({term_in_floor})")
    print(f"Bound = {m1_bound} + {m2_bound}")
    print(f"Bound = {total_bound}")
    print("\n---")
    print(f"The calculated upper bound is: {total_bound}")
    
    # Returning the value for the final answer tag
    return total_bound

# --- Parameters ---
# You can change these values to see the bound for different scenarios.
num_experts = 10
mistake_limit = 5

# Execute the calculation and print the result
final_bound = calculate_mistake_bound(num_experts, mistake_limit)
# The final answer tag will be based on the result of this execution.
# For n=10, c=5, the result is 26.
# For n=3, c=10, the result is 19.
# For n=101, c=20, the result is 1019.