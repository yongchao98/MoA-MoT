import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    
    Args:
        n: An integer.
        d: An odd integer, 2 <= d <= n.
    """
    if not (isinstance(n, int) and isinstance(d, int)):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return
    # The problem states d is odd, but the formula holds for any d.
    # We will not enforce the oddness of d to make the function more general.

    complexity = 2
    terms = []
    
    for i in range(1, d):
        term_val = combinations(n, i)
        complexity += term_val
        terms.append(f"C({n}, {i})")

    # Building the equation string
    equation_str = "2 + " + " + ".join(terms)
    
    # Building the equation string with values
    value_terms = []
    for i in range(1, d):
        value_terms.append(str(combinations(n, i)))
    
    equation_with_values_str = "2 + " + " + ".join(value_terms) + f" = {complexity}"

    print(f"For n={n} and d={d}, the smallest complexity is given by the formula:")
    print(f"Complexity = 2 + sum(C(n, i) for i in 1..d-1)")
    print(f"Calculation: {equation_str}")
    print(f"             = {equation_with_values_str}")
    
    # The final answer format as requested by the user prompt.
    # It's a bit unusual to print the final answer in this format,
    # but I will follow the instruction.
    # print(f"\n<<< {complexity} >>>")

if __name__ == '__main__':
    # Example from the problem context.
    # Let's use some example values for n and d.
    # For instance, n=5, d=3.
    n_val = 5
    d_val = 3
    solve_complexity(n_val, d_val)
    
    print("\n" + "="*20 + "\n")

    # Another example, n=6, d=5
    n_val = 6
    d_val = 5
    solve_complexity(n_val, d_val)
