import math

def calculate_min_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_{n,d}.

    Args:
        n (int): The upper range of indices for variables.
        d (int): The number of variables in each monomial.

    Returns:
        A tuple containing the final complexity and the list of term values.
    """
    if not (2 <= d <= n):
        raise ValueError("The constraints 2 <= d <= n must be met.")

    # The complexity is 2 + sum_{k=1}^{d-1} C(n, k)
    terms = []
    for k in range(1, d):
        try:
            term = math.comb(n, k)
            terms.append(term)
        except ValueError:
            # This can happen if k > n, but our check should prevent it.
            return -1, [] # Indicate an error

    complexity = 2 + sum(terms)
    return complexity, terms

def main():
    """
    Main function to demonstrate the calculation.
    """
    # Example values for n and d.
    # As per the problem, 2 <= d <= n and d is odd.
    n = 10
    d = 5 # d is odd

    print(f"Calculating the smallest complexity for n = {n} and d = {d}.")
    print(f"The formula for the complexity is C = 2 + sum_{k=1}^{d-1} C(n, k).")
    
    complexity, terms = calculate_min_complexity(n, d)

    if complexity != -1:
        # Build the equation string
        equation_str = "C = 2"
        for term in terms:
            equation_str += f" + {term}"
        
        # Add the final result
        equation_str += f" = {complexity}"

        print("The calculation is as follows:")
        print(equation_str)
    else:
        print("Could not compute complexity due to invalid inputs.")

if __name__ == "__main__":
    main()
