import sympy

def solve_limit_analytically():
    """
    This function calculates the specified limit using symbolic mathematics.
    It demonstrates the final step of the analytical derivation.
    """
    n = sympy.Symbol('n')
    
    # Based on the Central Limit Theorem, the probability P(n) is approximated by
    # P(n) ~= 1 - exp(-2/n) for large n.
    # We want to find the limit of n * P(n) as n -> infinity.
    
    one = 1
    minus_two = -2
    
    # Construct the expression for P(n)
    P_n = one - sympy.exp(minus_two/n)
    
    # The full expression for the limit
    expression = n * P_n
    
    # Calculate the limit using sympy
    limit_value = sympy.limit(expression, n, sympy.oo)
    
    print("The analytical derivation leads to the expression for the limit:")
    # Printing each numerical component of the equation as requested
    print(f"lim_{n}->oo  {n} * ({one} - exp({minus_two}/{n}))")
    
    print("\nThe result of this limit is:")
    print(limit_value)

if __name__ == '__main__':
    solve_limit_analytically()