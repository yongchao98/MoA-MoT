import sympy as sp

def solve_nonlinear_frequency_correction():
    """
    This function calculates the coefficient of the effective cubic nonlinearity, K_cubic,
    and extracts the coefficient of gamma^3 as the requested "3rd term".
    """
    # Define symbols
    gamma = sp.Symbol('gamma')
    A = sp.Symbol('A')
    A_bar = sp.Symbol('A_bar')
    
    # Linear oscillation frequency squared
    omega0_sq = 3 * gamma

    # Coefficients of the quadratic terms in the equation for x
    # from (3*gamma*(3*gamma+1)/2)*x**2 - x*x_ddot - (3/2)*x_dot**2
    # evaluated at the linear level (x_ddot = -omega0_sq * x)
    # The effective quadratic term is N1 = alpha * x0**2 + beta * x0_dot**2
    alpha = sp.S(3)/2 * gamma * (3*gamma + 1) + omega0_sq
    beta = -sp.S(3)/2

    # The solution to the O(epsilon) equation contains a constant term and a 2*omega0 term.
    # x1 = x1_0 + x1_2 * exp(2*i*omega0*t0) + c.c.
    # The complex amplitudes of these terms are:
    
    # C0 is the coefficient of the constant term in N1(x0)
    # x0**2 -> 2*A*A_bar, x0_dot**2 -> 2*omega0_sq*A*A_bar
    C0 = alpha * 2 * A * A_bar + beta * 2 * omega0_sq * A * A_bar
    x1_0 = C0 / omega0_sq
    
    # C2 is the coefficient of the exp(2*i*omega0*t0) term in N1(x0)
    # x0**2 -> A**2, x0_dot**2 -> -omega0_sq*A**2
    C2 = alpha * A**2 + beta * (-omega0_sq * A**2)
    x1_2 = C2 / (omega0_sq - (2*sp.sqrt(omega0_sq))**2)
    x1_2 = C2 / (omega0_sq - 4*omega0_sq)
    x1_2 = C2 / (-3 * omega0_sq)

    # The effective cubic nonlinearity K_cubic comes from several sources at O(epsilon^2):
    # 1. Interaction of quadratic terms: (2*alpha*x0*x1 + 2*beta*x0_dot*x1_dot)
    # 2. Explicit cubic term from R**(-3*gamma) expansion: C4 * x0**3
    
    # Contribution from x0*x1 term
    # The secular part is proportional to A*x1_0 + A_bar*x1_2
    term1_coeff = 2 * alpha * (A * x1_0 + A_bar * x1_2)
    
    # Contribution from x0_dot*x1_dot term
    # The secular part is proportional to (i*omega0*A)*(2*i*omega0*x1_2) + (-i*omega0*A_bar)*(-2*i*omega0*x1_2_bar)
    # which simplifies to 2 * beta * (2 * omega0_sq * A_bar * x1_2)
    term2_coeff = 2 * beta * (2 * omega0_sq * A_bar * x1_2)

    # Contribution from the explicit cubic term in the expansion of R**(-3*gamma)
    # The coefficient of epsilon^3*x^3 is -(3*gamma*(3*gamma+1)*(3*gamma+2))/6
    C4 = -(3*gamma * (3*gamma + 1) * (3*gamma + 2)) / 6
    # The secular part of x0**3 is 3 * A**2 * A_bar
    term3_coeff = C4 * 3 * A**2 * A_bar

    # Total coefficient K_cubic * A**2 * A_bar
    K_cubic_expr = term1_coeff.subs(A_bar, 1).subs(A, 1) + \
                   term2_coeff.subs(A_bar, 1).subs(A, 1) + \
                   term3_coeff.subs(A_bar, 1).subs(A, 1)
    
    # Simplify the expression for K_cubic
    K_cubic = sp.simplify(K_cubic_expr)
    
    # Expand the polynomial in gamma to find the coefficients
    K_cubic_poly = sp.Poly(K_cubic, gamma)
    
    # The "3rd term" is interpreted as the coefficient of gamma^3
    coeff_gamma_3 = K_cubic_poly.coeff_monomial(gamma**3)
    
    print(f"The effective cubic nonlinearity coefficient, K_cubic, is: {K_cubic}")
    print(f"The polynomial form of K_cubic is: {K_cubic_poly.as_expr()}")
    print(f"The coefficient of gamma^3 in K_cubic is the '3rd term'.")
    print(f"The 3rd term of the nonlinear correction is: {coeff_gamma_3}")

solve_nonlinear_frequency_correction()