import sympy

def solve_problem():
    """
    This script calculates the limit of the probability p_n based on a known result
    in probability theory for conditioned random walks.
    """
    
    # Define the distance and radius in terms of n
    d_n = "n"
    r_n = "n**(1/3)"
    
    # The formula for the limiting probability
    # lim_{n->inf} p_n = lim_{n->inf} [log(r(n)) / log(d(n))]
    
    # We will use sympy to represent the limit calculation
    n = sympy.Symbol('n')
    r = n**(sympy.Rational(1, 3))
    d = n
    
    expression = sympy.log(r) / sympy.log(d)
    
    # Calculate the limit
    limit_value = sympy.limit(expression, n, sympy.oo)
    
    # Print the steps of the calculation as requested
    print("The problem asks for the limit of p_n as n -> infinity.")
    print("p_n is the probability that a 2D random walk, conditioned to avoid the origin, hits a disk of radius r(n) centered at a distance d(n).")
    print("")
    print("According to a result by Benjamini, Lyons, and Schramm, this limit can be calculated with the formula:")
    print("lim_{n->inf} [log(r(n)) / log(d(n))]")
    print("")
    print(f"In this problem, the distance d(n) = {d_n} and the radius r(n) = {r_n}.")
    print("Substituting these into the formula, we get:")
    print(f"lim_{n->inf} [log({r_n}) / log({d_n})]")
    print("")
    print("Using the logarithm property log(x^a) = a*log(x), this simplifies to:")
    print("lim_{n->inf} [(1/3 * log(n)) / log(n)]")
    print("")
    print("The log(n) terms cancel out, leaving:")
    final_eq_val_1 = "1"
    final_eq_val_3 = "3"
    print(f"{final_eq_val_1}/{final_eq_val_3}")
    print("")
    
    # Final numeric result
    print(f"The calculated limit is: {float(limit_value)}")

solve_problem()
