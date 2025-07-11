import math

def get_mistake_bound(n, c):
    """
    Calculates the upper bound on the number of mistakes for the described
    experts problem variant.
    """

    # Helper function to calculate h(k) recursively.
    # h(k) represents the number of "generations" of majorities
    # that can be eliminated from a pool of k experts.
    memo = {}
    def h(k):
        if k in memo:
            return memo[k]
        if k <= 1:
            return 0
        
        # In each "generation", a minimal majority is eliminated. The size of
        # the remaining pool is ceil(k/2) - 1.
        next_k = math.ceil(k / 2.0) - 1
        result = 1 + h(next_k)
        memo[k] = result
        return result

    # 1. Decompose the total mistakes M into M1 and M2.
    #    M1: Algorithm and true expert are wrong.
    #    M2: Algorithm is wrong, but true expert is correct.

    # 2. Calculate the bound for M1.
    #    The true expert makes strictly fewer than c mistakes.
    #    So, the number of mistakes in M1 is at most c - 1.
    m1_bound = c - 1

    # 3. Calculate the bound for M2.
    #    There are n-1 "false" experts.
    n_false = n - 1
    
    # M2 is bounded by c * h(n-1).
    if n_false <= 0:
        h_val = 0
    else:
        h_val = h(n_false)
    
    m2_bound = c * h_val
    
    # 4. Combine the bounds for the total mistake bound.
    total_bound = m1_bound + m2_bound

    # Print the step-by-step explanation of the calculation.
    print(f"To find the upper bound on mistakes for n={n} experts and c={c} mistakes for elimination:")
    print("-" * 40)
    
    print("Step 1: Bound mistakes where the true expert is also wrong (M_1).")
    print(f"The true expert makes fewer than {c} mistakes.")
    print(f"Therefore, M_1 <= c - 1 = {m1_bound}.")
    print("-" * 40)

    print("Step 2: Bound mistakes where the true expert is correct but outvoted (M_2).")
    print(f"This involves the n - 1 = {n_false} false experts.")
    print("The bound is c * h(n-1), where h(k) is a recursive function representing rounds of eliminations.")
    print(f"For k = {n_false}, h({n_false}) = {h_val}.")
    print(f"Therefore, M_2 <= c * h(n-1) = {c} * {h_val} = {m2_bound}.")
    print("-" * 40)

    print("Step 3: Combine the bounds for the total.")
    print("Total Mistakes <= (Bound for M_1) + (Bound for M_2)")
    print(f"Total Mistakes <= ({c} - 1) + ({c} * {h_val})")
    print(f"Total Mistakes <= {m1_bound} + {m2_bound}")
    print(f"Total Mistakes <= {total_bound}")
    
    return total_bound

# --- User-configurable values ---
# Number of experts
n = 15
# Number of mistakes before an expert is removed
c = 10
# ------------------------------------

final_bound = get_mistake_bound(n, c)

print(f"\nThe final upper bound is {final_bound}.")

# The problem asks for the upper bound.
# The code calculates it and prints it as the final line.
# So we output this value in the requested format.
print(f"<<<{final_bound}>>>")