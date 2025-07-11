import sympy

def find_critical_exponent_order():
    """
    This function symbolically derives the order of the first correction
    to the critical exponent nu in phi^4 theory.
    """
    # Define the coupling constant u and a generic constant A as symbolic variables.
    u = sympy.Symbol('u')
    A = sympy.Symbol('A')

    # Step 1: State the relationship between the critical exponent nu and the
    # anomalous dimension of the mass term, gamma_m(u).
    # The mean-field value (for u=0) is nu = 1/2.
    print("The critical exponent nu is related to the coupling 'u' via the anomalous dimension gamma_m(u).")
    print("The formula is: 1 / nu = 2 + gamma_m(u)")
    print("For u=0, gamma_m(0)=0, giving the mean-field value nu = 1/2.")
    print("-" * 50)

    # Step 2: Define the series for gamma_m(u).
    # From one-loop calculations, the first contribution to gamma_m(u) is linear in 'u'.
    # We represent this as gamma_m(u) = A*u + O(u**2), where A is a non-zero constant.
    gamma_m_series = A * u + sympy.O(u**2)
    print(f"The anomalous dimension, from one-loop perturbation theory, expands as:")
    print(f"gamma_m(u) = {gamma_m_series}")
    print("-" * 50)

    # Step 3: Substitute gamma_m(u) into the expression for nu and expand.
    # We want to find the expansion of nu(u) = 1 / (2 + gamma_m(u)).
    nu_expr = 1 / (2 + A * u)
    
    # Perform a series expansion around u = 0 to find the correction terms.
    # The expansion is of nu(u) in powers of u.
    nu_series = nu_expr.series(u, 0, 3)

    print("Substituting this into the formula for nu and expanding in 'u' gives the final equation:")
    print(f"nu(u) = {nu_series}")
    print("-" * 50)

    # Step 4: Identify the order of the first correction term.
    # The expansion is nu(u) = 1/2 - (A/4)*u + O(u**2).
    # The term '-A*u/4' is the first non-vanishing contribution beyond the constant 1/2.
    # This term is proportional to u^1.
    print("The expansion shows the mean-field value of 1/2 and the subsequent correction terms.")
    print("The first correction term is (-A/4)*u.")
    print("This contribution is of first order in the coupling constant 'u'.")

find_critical_exponent_order()
