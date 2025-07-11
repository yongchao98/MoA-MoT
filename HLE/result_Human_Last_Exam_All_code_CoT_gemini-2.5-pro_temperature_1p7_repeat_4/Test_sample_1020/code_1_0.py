import math

def solve():
    """
    Calculates the smallest complexity for computing the polynomial f_n,d.
    """
    # Let n, d be integers with 2 <= d <= n and d odd.
    # We choose example values for n and d.
    n = 10
    d = 5

    # The complexity is given by the formula C = 2 + sum_{k=1}^{d-1} C(n, k)
    # where C(n,k) is the binomial coefficient "n choose k".
    
    # Let's check the constraints for our example values.
    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n is not met for n={n}, d={d}.")
        return
    if d % 2 == 0:
        print(f"Warning: The problem states d should be odd, but d={d} is even. The formula still applies.")

    # Calculate the sum of binomial coefficients
    terms = [math.comb(n, k) for k in range(1, d)]
    
    # Calculate the final complexity
    complexity = 2 + sum(terms)

    # Output the full equation as requested
    term_strings = [str(t) for t in terms]
    equation = f"2 + {' + '.join(term_strings)}"
    print(f"For n={n} and d={d}, the smallest complexity is given by the equation:")
    print(f"{equation} = {complexity}")

solve()