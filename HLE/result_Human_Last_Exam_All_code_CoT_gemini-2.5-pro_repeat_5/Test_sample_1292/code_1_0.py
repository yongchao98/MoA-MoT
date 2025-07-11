import sympy

def display_field_equation():
    """
    This function displays the derived field equation for the specified theory of gravity.
    Each term is printed separately to show the structure of the equation.
    """

    # Define the symbols for the tensors and scalars
    g, P, Q_tensor, Q_scalar, T = sympy.symbols(r'\sqrt{-g} P^{\alpha}_{\mu\nu} Q_{\nu}^{\alpha\beta} Q \mathrm{T}_{\mu\nu}')
    alpha, mu, nu = sympy.symbols(r'\alpha \mu \nu', cls=sympy.Function)

    # LHS terms
    term1 = sympy.sympify("-2/sqrt(-g) * d/d_alpha * (sqrt(-g) * P_alpha_mu_nu)", 
                          locals={'P_alpha_mu_nu': P, 'sqrt': sympy.sqrt, 'g': sympy.Symbol('-g'), 'd_alpha': sympy.Symbol(r'\partial_{\alpha}')})
    term2 = sympy.sympify("-2 * P_mu_alpha_beta * Q_nu_alpha_beta",
                          locals={'P_mu_alpha_beta': sympy.Symbol(r'P_{\mu\alpha\beta}'), 'Q_nu_alpha_beta': sympy.Symbol(r'Q_{\nu}^{\alpha\beta}')})
    term3 = sympy.sympify("Q_alpha_beta_mu * P_alpha_beta_nu",
                          locals={'Q_alpha_beta_mu': sympy.Symbol(r'Q^{\alpha\beta}_{\mu}'), 'P_alpha_beta_nu': sympy.Symbol(r'P_{\alpha\beta\nu}')})
    term4 = sympy.sympify("-1/2 * Q_scalar * g_mu_nu",
                          locals={'Q_scalar': Q_scalar, 'g_mu_nu': sympy.Symbol(r'g_{\mu\nu}')})

    # RHS term
    rhs = sympy.sympify("(8*pi*G/c**4) * T_mu_nu",
                        locals={'pi': sympy.pi, 'G': sympy.Symbol('G'), 'c': sympy.Symbol('c'), 'T_mu_nu': T})

    # Print the equation part by part
    print("The derived field equation is:")
    
    # Using string formatting for better alignment and explicit representation
    term1_str = f"-\u2202\u03b1(\u221A-g P\u03b1\u03bc\u03bd) / \u221A-g"
    term2_str = f"- 2 P\u03bc\u03b1\u03b2 Q\u03bd\u03b1\u03b2"
    term3_str = f"+ Q\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd"
    term4_str = f"- (1/2) Q g\u03bc\u03bd"
    rhs_str = f"= (8\u03c0G/c\u2074) T\u03bc\u03bd"

    # Simplified representation for clarity. This does not fully capture the tensor indices in SymPy symbols.
    # The string representation is more direct for this kind of output.
    print(f"Term 1 (Derivative of Superpotential): -2/\u221A-g \u2202\u03b1(\u221A-g P\u03b1\u03bc\u03bd)")
    print(f"Term 2: - 2 P\u03bc\u03b1\u03b2 Q\u03bd\u03b1\u03b2")
    print(f"Term 3: + Q\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd")
    print(f"Term 4 (Non-metricity scalar term): - (1/2) Q g\u03bc\u03bd")
    print(f"Equals: =")
    print(f"RHS (Energy-Momentum Tensor): (8\u03c0G/c\u2074) T\u03bc\u03bd")
    
    print("\nFull Equation:")
    print(f" -2/\u221A-g \u2202\u03b1(\u221A-g P\u03b1\u03bc\u03bd) - 2 P\u03bc\u03b1\u03b2 Q\u03bd\u03b1\u03b2 + Q\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd - (1/2) Q g\u03bc\u03bd = (8\u03c0G/c\u2074) T\u03bc\u03bd")


display_field_equation()