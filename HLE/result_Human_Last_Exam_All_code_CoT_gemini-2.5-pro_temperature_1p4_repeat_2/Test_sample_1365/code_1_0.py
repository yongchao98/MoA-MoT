import math

def calculate_mistake_bound(n, c):
    """
    Calculates and prints the upper bound on the number of mistakes for the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The mistake threshold for removing an expert.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be positive integers.")
        return

    # Bound for M'': mistakes where the true expert is also wrong.
    # The true expert makes at most c-1 mistakes.
    m_double_prime_bound = c - 1
    
    # Bound for M': mistakes where the true expert is right.
    # As derived, M' <= (n-1)*c / 2.
    m_prime_bound_float = (n - 1) * c / 2
    m_prime_bound = math.floor(m_prime_bound_float)

    # The total bound is the sum of the two cases.
    total_bound = m_prime_bound + m_double_prime_bound

    print(f"For n = {n} experts and a mistake tolerance c = {c}:")
    print("-" * 40)
    print("The total mistake bound is the sum of two cases:")
    print("1. M'' (true expert is wrong):")
    print(f"   Bound = c - 1 = {c} - 1 = {m_double_prime_bound}")
    print("\n2. M' (true expert is right):")
    print(f"   Bound = floor((n-1) * c / 2) = floor(({n}-1) * {c} / 2) = floor({m_prime_bound_float}) = {m_prime_bound}")
    print("-" * 40)
    print(f"Total Upper Bound M <= M' + M''")
    print(f"Total Upper Bound M <= {m_prime_bound} + {m_double_prime_bound} = {total_bound}")
    
# Example: Provide values for n and c to calculate the bound.
# You can change these values to see the bound for different scenarios.
n_experts = 21
c_tolerance = 10

calculate_mistake_bound(n_experts, c_tolerance)

# Let's calculate the bound for the user's specific (but unstated) n and c
# Since they are not given, let's assume they might be, e.g. n=5, c=3
# This is just an example calculation for a different pair of values.
# To run your own, modify n_experts and c_tolerance above.
# print("\n--- Another Example ---")
# calculate_mistake_bound(5, 3)
