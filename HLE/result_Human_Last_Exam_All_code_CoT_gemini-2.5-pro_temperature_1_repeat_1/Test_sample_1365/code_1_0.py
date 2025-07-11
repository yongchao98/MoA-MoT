import math

def calculate_mistake_bound(n, c):
    """
    Calculates an upper bound on the number of mistakes for a majority voting algorithm.

    Args:
        n (int): The total number of experts.
        c (int): The number of mistakes an expert must make to be removed.
    """
    if not isinstance(n, int) or not isinstance(c, int) or n <= 1 or c < 1:
        print("Please provide positive integers for n (n>1) and c (c>=1).")
        return

    # M_2 is the bound on mistakes when the true expert is wrong.
    # The true expert makes strictly fewer than c mistakes.
    m2_bound = c - 1

    # M_1 is the bound on mistakes when the true expert is correct.
    # This involves a logarithmic number of phases, each taking c mistakes.
    if n <= 2:
        # If there's only one non-true expert, M_1 is 0 as there's no majority among them.
        m1_bound_log_term = 0
    else:
        # The base of the logarithm is 1.5, from the pool reduction factor.
        m1_bound_log_term = math.log(n - 1) / math.log(1.5)
    
    m1_bound = c * m1_bound_log_term
    
    # Total bound is the sum of the two bounds.
    total_bound = m1_bound + m2_bound

    print("Derivation of the upper bound:")
    print(f"n = {n}")
    print(f"c = {c}")
    print("-" * 30)
    
    print(f"1. Bound on mistakes when true expert is wrong (M_2):")
    print(f"   M_2 <= c - 1 = {c} - 1 = {m2_bound}")
    print("-" * 30)

    print(f"2. Bound on mistakes when true expert is correct (M_1):")
    print(f"   M_1 <= c * log_1.5(n-1)")
    print(f"   M_1 <= {c} * log_1.5({n-1})")
    print(f"   M_1 <= {c} * {m1_bound_log_term:.4f}")
    print(f"   M_1 <= {m1_bound:.4f}")
    print("-" * 30)

    print("3. Total Upper Bound (M <= M_1 + M_2):")
    # Using the calculated values in the final print statement
    final_equation_str = f"M <= ({c} * {m1_bound_log_term:.4f}) + {m2_bound}"
    print(f"   {final_equation_str}")
    print(f"   M <= {m1_bound:.4f} + {m2_bound}")
    print(f"   M <= {total_bound:.4f}")

# Example usage with n=10 experts and c=5 mistakes threshold.
# Feel free to change these values.
n_experts = 10
c_mistakes = 5
calculate_mistake_bound(n_experts, c_mistakes)