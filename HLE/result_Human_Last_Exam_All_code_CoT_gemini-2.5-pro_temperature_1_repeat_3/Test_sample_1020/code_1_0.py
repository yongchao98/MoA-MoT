import math

def solve():
    """
    Calculates the smallest complexity of a matrix product for f_n,d.
    """
    # Set the integer parameters n and d.
    # n must be greater than or equal to d.
    # d must be odd and at least 2.
    n = 10
    d = 5

    # --- Parameter validation ---
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if d < 2 or d > n:
        print(f"Error: Constraints not met. Required: 2 <= d <= n. Got n={n}, d={d}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd. Got d={d}.")
        return

    # --- Complexity calculation ---
    # The minimal complexity is C = 2 + sum_{i=1}^{d-1} max(C(n,i), C(n,d-i))
    
    terms_values = []
    for i in range(1, d):
        # Calculate the binomial coefficients C(n, i) and C(n, d-i)
        try:
            comb_i = math.comb(n, i)
            comb_di = math.comb(n, d - i)
        except ValueError:
            print(f"Error: math.comb(n, k) failed. n={n}, k is {i} or {d-i}. k must be between 0 and n.")
            return

        # The size of the i-th cut is the maximum of these two
        term = max(comb_i, comb_di)
        terms_values.append(term)
        
    total_complexity = 2 + sum(terms_values)
    
    # --- Output the equation as requested ---
    # Example: C = 2 + 6 + 6 = 14
    
    # Create the string for the sum part, e.g., "6 + 6"
    sum_str = " + ".join(map(str, terms_values))
    
    # Create the final equation string
    equation = f"C = 2 + {sum_str} = {total_complexity}"
    
    print(f"For n={n} and d={d}, the smallest complexity is calculated as:")
    print(equation)

solve()