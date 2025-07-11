import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    The user will be prompted to enter the integers n and d.
    """
    try:
        n_str = input("Enter the integer n: ")
        d_str = input("Enter the integer d: ")
        n = int(n_str)
        d = int(d_str)

        if not (2 <= d <= n):
            raise ValueError("Constraints not met: must have 2 <= d <= n.")
        if d % 2 == 0:
            raise ValueError("Constraint not met: d must be an odd integer.")

    except ValueError as e:
        print(f"Invalid input: {e}")
        return

    # m = (d-1)/2
    m = (d - 1) // 2

    sum_terms_str = []
    sum_terms_val = []

    for k in range(1, m + 1):
        # Calculate binomial coefficients C(n,k) and C(n,d-k)
        try:
            c1 = math.comb(n, k)
            c2 = math.comb(n, d - k)
        except ValueError:
            # This can happen if k or d-k are > n, though our loops prevent k>n.
            # If d>n, we would have already exited.
            # So this is just for safety.
            c1 = 0
            c2 = 0

        # Add strings for formula display
        sum_terms_str.append(f"min(C({n},{k}), C({n},{d-k}))")
        
        # Add numeric values for calculation
        sum_terms_val.append(min(c1, c2))

    # Calculate final complexity
    rank_sum = sum(sum_terms_val)
    complexity = 2 + 2 * rank_sum

    # Print the step-by-step calculation
    print("\nThe complexity C is given by the formula:")
    print(f"C = 2 + 2 * sum_{{k=1}}^{{{m}}}[min(C(n,k), C(n,d-k))]")
    print("\nPlugging in n =", n, "and d =", d, ":")
    
    if not sum_terms_str:
        # This case happens if m=0, i.e., d=1, but we have d>=2.
        # Or if d=3, m=1. Let's handle all cases gracefully.
        # For d=3, m=1, loop runs once.
        # So sum_terms_str will not be empty if d>=3.
        # In case m=0 (d=1), the sum is empty (0).
        formula_str = "0"
        values_str = "0"
    else:
        formula_str = " + ".join(sum_terms_str)
        values_str = " + ".join(map(str, sum_terms_val))

    print(f"C = 2 + 2 * ({formula_str})")
    print(f"C = 2 + 2 * ({values_str})")
    if len(sum_terms_val) > 1:
        print(f"C = 2 + 2 * ({rank_sum})")
    print(f"C = 2 + {2 * rank_sum}")
    print(f"C = {complexity}")
    
    # Returning the final answer in the requested format
    print("\nFinal Answer:")
    print(f"<<<{complexity}>>>")


if __name__ == "__main__":
    solve_complexity()