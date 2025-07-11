def calculate_mistake_bound(n, c):
    """
    Calculates and explains the upper bound on mistakes for a majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 0 or c <= 0:
        print("Error: n and c must be positive integers.")
        return

    # Upper bound on mistakes when the true expert is wrong
    # This is at most c-1, as the true expert makes < c mistakes.
    bound_m_incorrect = c - 1

    # Upper bound on mistakes when the true expert is correct
    # This is at most the total number of mistakes the n-1 false experts can make.
    bound_m_correct = (n - 1) * c

    # The total upper bound is the sum of the two bounds.
    total_bound = bound_m_correct + bound_m_incorrect

    # The simplified formula for the total bound is n*c - 1.
    simplified_bound = n * c - 1

    print(f"For n = {n} experts and a mistake threshold c = {c}:")
    print("\nThe upper bound on algorithm mistakes is the sum of two cases:")
    
    # Print the equation for M_incorrect
    print("\n1. Max mistakes when the true expert is wrong (M_incorrect):")
    print(f"   M_incorrect <= c - 1")
    print(f"   M_incorrect <= {c} - 1 = {bound_m_incorrect}")

    # Print the equation for M_correct
    print("\n2. Max mistakes when the true expert is correct (M_correct):")
    print(f"   M_correct <= (n - 1) * c")
    print(f"   M_correct <= ({n} - 1) * {c} = {n - 1} * {c} = {bound_m_correct}")

    # Print the equation for the total bound
    print("\nTotal Upper Bound = M_incorrect + M_correct")
    print(f"Total Upper Bound <= {bound_m_incorrect} + {bound_m_correct} = {total_bound}")
    
    # Verify with the simplified formula
    print("\nThis simplifies to the formula: n * c - 1")
    print(f"n * c - 1 = {n} * {c} - 1 = {n*c} - 1 = {simplified_bound}")

    # The two calculation methods should yield the same result
    assert total_bound == simplified_bound

# --- Example Usage ---
# You can change these values to see the bound for different scenarios.
num_experts = 20
mistake_threshold = 10
calculate_mistake_bound(num_experts, mistake_threshold)