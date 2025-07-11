import math

def calculate_smallest_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing the polynomial f_n,d.

    Args:
        n (int): The upper limit for the indices of the variables.
        d (int): The degree of the monomials in the polynomial.
    """
    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n must be met. Received n={n}, d={d}.")
        return

    # The formula for the smallest complexity is 2 + sum(C(n, i) for i from 1 to d-1)
    complexity = 2
    terms = []

    for i in range(1, d):
        try:
            term = math.comb(n, i)
            terms.append(term)
            complexity += term
        except ValueError:
            print(f"Error: Cannot compute C({n}, {i}). Check if n >= i.")
            return

    # Building the explanation string
    formula_str = "2 + " + " + ".join([f"C({n}, {i+1})" for i in range(len(terms))])
    values_str = "2 + " + " + ".join(map(str, terms))

    print(f"For n={n} and d={d}:")
    print(f"The smallest complexity is given by the formula: {formula_str}")
    print(f"Calculation: {values_str} = {complexity}")
    
    return complexity


if __name__ == '__main__':
    # You can change these values to compute the complexity for different n and d.
    # Example values as per the problem context:
    # Let's use n=10, d=5 (d is odd and 2 <= d <= n)
    n_example = 10
    d_example = 5
    
    calculate_smallest_complexity(n_example, d_example)

    print("\n" + "="*30 + "\n")

    # Another example
    n_example_2 = 5
    d_example_2 = 3
    calculate_smallest_complexity(n_example_2, d_example_2)
