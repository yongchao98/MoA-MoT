import math

def calculate_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for a majority voting algorithm
    in a variant of the experts problem.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert can make before being removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n < 1 or c < 1:
        print("Error: n and c must be integers greater than or equal to 1.")
        return

    # The total number of mistakes, M, is the sum of two types of mistakes:
    # M_incorr: Algorithm mistakes where the true expert is also incorrect.
    # M_corr: Algorithm mistakes where the true expert is correct.
    # M = M_incorr + M_corr

    # 1. Bound for M_incorr
    # The true expert makes strictly fewer than c mistakes. M_incorr is a subset
    # of these mistakes, so M_incorr has an upper bound of c - 1.
    bound_m_incorr = c - 1

    # 2. Bound for M_corr
    # For each M_corr mistake, at least 2 non-true experts must have been wrong.
    # The total mistake "budget" for the (n-1) non-true experts is (n-1)*c.
    # So, 2 * M_corr <= (n-1)*c  =>  M_corr <= (n-1)*c / 2
    # Since M_corr must be an integer, we take the floor.
    bound_m_corr = math.floor((n - 1) * c / 2)

    # 3. Combine the bounds
    total_bound = bound_m_corr + bound_m_incorr

    print(f"Calculating the upper bound for n={n} experts and mistake threshold c={c}:")
    print("-" * 20)
    print("The final formula for the upper bound M is: floor((n - 1) * c / 2) + (c - 1)")
    print("\nStep-by-step calculation:")
    
    # Output each number in the final equation
    print(f"Bound for mistakes where true expert is incorrect (M_incorr) <= c - 1")
    print(f"   = {c} - 1 = {bound_m_incorr}")

    print(f"\nBound for mistakes where true expert is correct (M_corr) <= floor((n - 1) * c / 2)")
    print(f"   = floor(({n} - 1) * {c} / 2)")
    print(f"   = floor({(n - 1)} * {c} / 2)")
    print(f"   = floor({(n - 1) * c} / 2)")
    print(f"   = floor({(n - 1) * c / 2}) = {bound_m_corr}")
    
    print(f"\nTotal Upper Bound M <= M_corr + M_incorr")
    print(f"   = {bound_m_corr} + {bound_m_incorr} = {total_bound}")
    print("-" * 20)
    print(f"The upper bound on the number of mistakes is: {total_bound}")


if __name__ == '__main__':
    # Example values for n and c.
    # You can change these to test other scenarios.
    num_experts = 11
    mistake_threshold = 5
    calculate_mistake_bound(num_experts, mistake_threshold)