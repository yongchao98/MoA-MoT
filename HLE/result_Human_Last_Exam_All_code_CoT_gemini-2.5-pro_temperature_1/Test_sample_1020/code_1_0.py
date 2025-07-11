import math

def solve():
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_n,d.
    """
    # Let n, d be integers with 2 <= d <= n and d odd.
    # We choose a specific example for n and d that satisfies the constraints.
    n = 5
    d = 3

    if not (2 <= d <= n and d % 2 != 0):
        print("The chosen n and d do not satisfy the problem constraints.")
        return

    # The smallest complexity is given by the formula: 2 + sum(C(n, i) for i in 1..d-1)
    # where C(n, i) is the binomial coefficient "n choose i".
    # We calculate the sum of the binomial coefficients.
    
    binomial_terms = []
    for i in range(1, d):
        try:
            term = math.comb(n, i)
            binomial_terms.append(term)
        except ValueError:
            print(f"Cannot compute C({n}, {i}). Check that n >= i.")
            return

    total_complexity = 2 + sum(binomial_terms)

    # Output the explanation and the final equation with all numbers.
    print(f"To find the smallest complexity for n={n} and d={d}, we use the formula:")
    print("Complexity = 2 + sum(C(n, i) for i from 1 to d-1)")
    print("\nThis is derived from a construction where matrix dimensions correspond to the number of ways to choose subsets of indices.")
    
    # Build the string for the equation to show the calculation.
    equation_str = "2 + " + " + ".join(map(str, binomial_terms))
    
    print(f"\nFor n={n}, d={d}, the calculation is:")
    print(f"Complexity = {equation_str}")
    print(f"Complexity = {total_complexity}")

solve()