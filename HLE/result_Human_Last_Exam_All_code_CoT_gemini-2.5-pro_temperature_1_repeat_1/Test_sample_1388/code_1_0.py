import sympy

def solve_pde_estimate():
    """
    This function determines the explicit form of H(t) in an L2 energy estimate
    for a given non-local PDE. It uses sympy to construct and print the formula.
    """
    
    # Define the symbols for the quantities in the formula
    h_t = sympy.Symbol('h(t)', real=True, positive=True)
    L1_norm_u0 = sympy.Symbol('||u_0||_L1', real=True, positive=True)
    L2_norm_u0_sq = sympy.Symbol('||u_0||_L2**2', real=True, positive=True)
    
    # Construct the expression for H(t)^2 based on the derivation
    # H(t)^2 = exp(-2*h(t)) + (||u_0||_L1 / ||u_0||_L2**2) * (1 - exp(-2*h(t)))
    
    # Each number in the equation is explicitly written out
    term1_exponent = -2 * h_t
    term1 = sympy.exp(term1_exponent)
    
    term2_factor = 1 - sympy.exp(-2 * h_t)
    term2 = (L1_norm_u0 / L2_norm_u0_sq) * term2_factor
    
    H_t_squared = term1 + term2
    
    # H(t) is the square root of the expression
    H_t = sympy.sqrt(H_t_squared)
    
    # Print the final expression for H(t)
    print("The explicit form of H(t) is:")
    sympy.pprint(H_t, use_unicode=False)

solve_pde_estimate()