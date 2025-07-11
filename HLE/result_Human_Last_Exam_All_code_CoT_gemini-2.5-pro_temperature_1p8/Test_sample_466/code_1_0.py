import math

def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category
    from F_11^3 to itself. This number is equal to the order of the
    general linear group GL(3, 11).
    """

    n = 3
    q = 11

    # The order of GL(n, q) is given by the product (q^n - q^0) * ... * (q^n - q^(n-1)).
    print(f"The number of adjunctions is the order of the general linear group GL({n}, {q}).")
    print("The formula for the order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.")
    print(f"For n={n} and q={q}, this is:")

    # Calculate q^n
    q_to_n = q**n

    # Print the formula with values substituted
    term_formulas = []
    for i in range(n):
        term_formulas.append(f"({q}^{n} - {q}^{i})")
    print(" * ".join(term_formulas))

    # Print the formula with q^n computed
    computed_term_formulas = []
    for i in range(n):
        computed_term_formulas.append(f"({q_to_n} - {q**i})")
    print("= " + " * ".join(computed_term_formulas))
    
    # Calculate each term and the final result
    terms = []
    result = 1
    for i in range(n):
        term = q_to_n - q**i
        terms.append(str(term))
        result *= term

    # Print the expanded equation
    print(f"= {' * '.join(terms)}")
    
    # Print the final result
    print(f"= {result}")

solve()