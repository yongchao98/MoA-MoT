import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    # Using math.comb is efficient and safe for Python 3.8+
    # For compatibility, a factorial-based implementation is also fine.
    try:
        return math.comb(n, k)
    except AttributeError:
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
    Calculates the smallest known complexity of a matrix product for f_n,d.
    """
    # You can change these values for n and d.
    # As per the problem description, d must be odd and 2 <= d <= n.
    n = 10
    d = 5

    print(f"Calculating the complexity for n = {n} and d = {d}:")

    if not (isinstance(n, int) and isinstance(d, int) and 2 <= d <= n and d % 2 != 0):
        print("Warning: The provided n and d values do not satisfy the constraints (2 <= d <= n, d is odd).")
        print("The calculation will proceed, but the theoretical grounding may not apply.")

    # The complexity is given by the formula: 2 + sum_{i=1}^{d-1} C(n, i)
    sum_of_combinations = 0
    
    # Building the string for the equation with symbolic C(n,k)
    equation_symbolic_parts = []
    for i in range(1, d):
        equation_symbolic_parts.append(f"C({n},{i})")
    equation_symbolic_str = "2 + " + " + ".join(equation_symbolic_parts)

    # Building the string for the equation with numerical values and calculating the sum
    equation_numeric_parts = []
    for i in range(1, d):
        comb_val = combinations(n, i)
        sum_of_combinations += comb_val
        equation_numeric_parts.append(str(comb_val))
    equation_numeric_str = "2 + " + " + ".join(equation_numeric_parts)

    complexity = 2 + sum_of_combinations

    print("\nThe smallest known complexity is given by the formula:")
    print(f"Complexity = {equation_symbolic_str}")
    print("\nSubstituting the values:")
    print(f"Complexity = {equation_numeric_str}")
    print(f"\nFinal calculated complexity:")
    print(f"Complexity = {complexity}")

solve_complexity()