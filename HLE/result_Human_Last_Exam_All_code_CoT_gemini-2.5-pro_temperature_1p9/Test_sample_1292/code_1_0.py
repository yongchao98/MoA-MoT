import sympy as sp

def print_field_equation():
    """
    This function prints the derived field equation for Symmetric Teleparallel Gravity.
    The equation is represented using sympy for clear mathematical formatting.
    """
    # Define indices
    mu, nu, alpha, beta, rho, sigma = sp.symbols('mu nu alpha beta rho sigma', cls=sp.Idx)
    
    # Define tensors and scalars
    g = sp.Function('g')      # Metric tensor
    P = sp.Function('P')      # Superpotential tensor
    Q = sp.Function('Q')      # Non-metricity tensor
    Q_scalar = sp.Symbol('Q') # Non-metricity scalar
    T = sp.Function('T')      # Energy-momentum tensor
    G = sp.Symbol('G')        # Gravitational constant
    c = sp.Symbol('c')        # Speed of light
    
    # Square root of negative determinant of the metric
    sqrt_neg_g = sp.sqrt(-sp.Symbol('g_det'))
    
    # First term: derivative of the superpotential
    term1 = -2 / sqrt_neg_g * sp.Derivative(sqrt_neg_g * P(alpha, mu, nu), alpha)
    
    # Second term: contraction of P and Q
    term2 = -2 * P(mu, alpha, beta) * Q(nu, alpha, beta, "up")
    
    # Third term: contraction of Q and P
    term3 = Q(alpha, beta, mu, "up") * P(alpha, beta, nu)
    
    # Fourth term: non-metricity scalar term
    term4 = -sp.Rational(1, 2) * Q_scalar * g(mu, nu)
    
    # Right-hand side: energy-momentum tensor
    rhs = (8 * sp.pi * G / c**4) * T(mu, nu)

    # Define the equation components for printing
    eq_parts = {
        "Term 1 (Derivative of P)": term1,
        "Term 2 (P*Q contraction)": term2,
        "Term 3 (Q*P contraction)": term3,
        "Term 4 (Scalar Q term)": term4,
        "RHS (Energy-Momentum)": rhs
    }

    print("The derived field equation is:")
    
    # We build the equation string manually for better formatting control,
    # as direct sympy printing can be hard to read for complex tensor equations.
    
    # Print each term in the final equation
    print(f"- (2/sqrt(-g)) * d/dx^alpha (sqrt(-g) * P^{{alpha}}_{{mu nu}})")
    print(f"- 2 * P_{{mu alpha beta}} * Q_{{nu}}^{{alpha beta}}")
    print(f"+ Q^{{alpha beta}}_{{mu}} * P_{{alpha beta nu}}")
    print(f"- (1/2) * Q * g_{{mu nu}}")
    print(f"= (8*pi*G/c^4) * T_{{mu nu}}")
    
    print("\nSymbolic representation of each number in the final equation:")
    print(f"First term coefficient: -2")
    print(f"Second term coefficient: -2")
    print(f"Third term coefficient: +1")
    print(f"Fourth term coefficient: -1/2")
    print(f"Right-hand side coefficient: 8")


print_field_equation()