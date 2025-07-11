import sympy

def derive_effective_interaction():
    """
    Performs symbolic derivation of the effective electron-electron interaction
    mediated by phonons in the static limit.
    """
    # Define symbolic variables
    g, q_j, m = sympy.symbols('g q_j m')
    w_q = sympy.Symbol('w_q', positive=True) # Phonon frequency
    rho_q = sympy.Symbol('rho_q')           # Electron density operator for mode q
    rho_minus_q = sympy.Symbol('rho_-q')   # Electron density operator for mode -q
    
    # 1. Define the coupling constant C_qj
    # We use its squared magnitude |C_qj|^2 as derived in the explanation
    C_qj_squared = (g**2 * q_j**2) / (2 * m * w_q)
    
    # 2. Define the static effective potential V_eff(q,j) = -2 * |C_qj|^2 / w_q
    V_eff = -2 * C_qj_squared / w_q
    
    # 3. Construct the effective interaction Hamiltonian for a single mode (q,j)
    # H_eff = 1/2 * V_eff * rho_q * rho_-q
    H_eff = sympy.Rational(1, 2) * V_eff * rho_q * rho_minus_q
    
    # Simplify the final expression
    H_eff_simplified = sympy.simplify(H_eff)
    
    # Print the result in a formatted equation
    print("The effective electron-electron interaction for a single mode (q, j) is:")
    
    # Create the equation string for printing
    lhs = "H_eff(q, j)"
    
    # Pretty print the final Hamiltonian term by term
    # Decompose the expression to print parts separately
    coeff, operators = H_eff_simplified.as_coeff_mul()
    
    # Get the numerator and denominator of the coefficient
    num, den = coeff.as_numer_denom()

    # Reconstruct the string for printing to match the equation format
    # The negative sign is pulled out front
    if num < 0:
        sign = "-"
        num = -num
    else:
        sign = ""
        
    term_str = f"{g**2 * q_j**2}"
    denom_str = f"{den}"
    
    # Assemble final output
    print(f"{lhs} = {sign} ({term_str}) / ({denom_str}) * {rho_q} * {rho_minus_q}")
    
    
derive_effective_interaction()
