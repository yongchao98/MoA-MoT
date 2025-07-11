import sympy as sp

def solve_susceptibility():
    """
    This function symbolically derives the magnetic susceptibility chi
    for an Ising model on a sparse random graph based on the provided formula and assumptions.
    """
    # Step 1: Define all symbols
    chi, beta, J, c, m0, N = sp.symbols('chi beta J c m_0 N')
    l = sp.Symbol('l', positive=True, integer=True)

    # Step 2: Derive the correlation function C_l
    # Based on the model m_i = tanh(beta*(B_i + J*sum(m_j))), the linear response analysis shows:
    # d<sigma_0>/dB_l = beta * (1 - m0**2) * (beta * J * (1 - m0**2))**l
    # C_l = (1/beta) * d<sigma_0>/dB_l
    
    # Let X = beta * J * (1 - m0**2) for simplicity
    X = beta * J * (1 - m0**2)
    # The correlation function C_l is then:
    C_l = (1 - m0**2) * X**l

    # Step 3: Compute the susceptibility chi by summing the series
    # The formula is chi = beta * Sum_{l=1 to inf} c*(c-1)**(l-1) * C_l
    # We can pull constants out of the sum:
    # chi = beta * c * (1-m0**2) * X * Sum_{l=1 to inf} [(c-1)*X]**(l-1)
    # This is a geometric series Sum_{k=0 to inf} r**k = 1/(1-r)
    # where r = (c-1)*X
    
    # So the sum evaluates to 1 / (1 - (c-1)*X)
    geometric_sum = 1 / (1 - (c - 1) * X)
    
    # Substitute back into the expression for chi
    chi_expr = beta * c * (1 - m0**2) * X * geometric_sum

    # Substitute X back to get chi in terms of basic variables
    chi_expr = chi_expr.subs(X, beta * J * (1 - m0**2))
    
    # Simplify the expression
    chi_expr = sp.simplify(chi_expr)
    
    # The intermediate result is: chi = (beta**2*c*J*(1 - m0**2)**2) / (1 - (c - 1)*beta*J*(1 - m0**2))

    # Step 4: Simplify using the definition of N
    # N = beta * c * (1 - m0**2) / (c - 1)
    # We solve for (1 - m0**2) to substitute it into the expression for chi
    
    eq_N = sp.Eq(N, beta * c * (1 - m0**2) / (c - 1))
    one_minus_m0_sq_expr = sp.solve(eq_N, (1 - m0**2))[0]
    
    # Substitute this into the expression for chi
    final_chi_expr = chi_expr.subs((1 - m0**2), one_minus_m0_sq_expr)
    
    # Simplify the final expression
    final_chi_expr = sp.simplify(final_chi_expr)

    # Print the final equation
    # The problem requires printing each "number" in the equation.
    # We will format the output to clearly show the final equation chi = ...
    print("The final expression for the magnetic susceptibility is:")
    
    # For clarity, let's build the equation step-by-step for printing
    numerator_terms = []
    denominator_terms = []

    # Deconstruct the final expression for printing
    num, den = final_chi_expr.as_numer_denom()
    
    num_str = str(num).replace('**','^')
    den_str = str(den).replace('**','^')
    
    print(f"\u03C7 = ({num_str}) / ({den_str})")


solve_susceptibility()

# Final Answer Wrapper
# The final derived formula is chi = (J * N**2 * (c - 1)**2) / (c - J * N * (c - 1)**2)
final_answer = 'J*N^2*(c - 1)^2 / (c - J*N*(c - 1)^2)'
print(f'<<<{final_answer}>>>')
