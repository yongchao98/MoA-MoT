def print_minimum_cost_formula():
    """
    This function prints the derived formula for the minimum achievable curvature cost.
    The cost is expressed in Big-O notation as a function of 'n' (number of samples)
    and 'd' (dimension of the space), where n < d.
    """

    # The derived minimum cost is the complexity of the most efficient algorithm
    # for the NGD update, which exploits the low-rank structure of the Fisher matrix.
    
    # The prompt requires printing the numbers in the final equation.
    # In the complexity formula O(n*d^2 + n^2*d), the numbers are the exponents.
    # We will print the formula to reflect this.

    n = 'n'
    d = 'd'
    power_2 = '2'

    print("The minimum achievable curvature cost is determined by the computational complexity of an efficient NGD update algorithm.")
    print("The complexity, in Big-O notation, is:")

    # Print the expanded form of the complexity: O(n*d^2 + n^2*d)
    # This clearly shows the constituent costs and their exponents.
    print(f"O({n}*{d}^{power_2} + {n}^{power_2}*{d})")

    # For completeness, also print the factored form.
    print("\nThis can be factored as:")
    print(f"O({n}*{d}*({n} + {d}))")

print_minimum_cost_formula()
