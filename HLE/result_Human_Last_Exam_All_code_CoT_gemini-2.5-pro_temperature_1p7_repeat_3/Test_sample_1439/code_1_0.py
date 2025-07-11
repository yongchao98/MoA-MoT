import sympy as sp

def analyze_critical_exponent_nu():
    """
    Analyzes the perturbative expansion of the critical exponent nu.

    Within the framework of the renormalization group for phi^4 theory, the critical
    exponent nu is related to the anomalous dimension of the phi^2 operator, gamma_phi2(u).

    This function symbolically expands nu(u) = 1 / (2 + gamma_phi2(u)) to show
    at which order in the coupling 'u' the first correction to the mean-field
    value (nu=1/2) appears.
    """
    # Define symbols for the coupling 'u' and the first expansion
    # coefficient 'c1' of the anomalous dimension.
    u = sp.Symbol('u')
    c1 = sp.Symbol('c1')

    # The anomalous dimension gamma_phi2(u) has a perturbative expansion in 'u'.
    # The first, one-loop contribution is linear in 'u'. We only need this
    # leading term to find the first correction to nu.
    # gamma_phi2(u) = c1*u + O(u**2)
    gamma_phi2 = c1 * u

    # The formula for the critical exponent nu is given by:
    nu = 1 / (2 + gamma_phi2)

    # Perform a Taylor series expansion of nu around u=0 to find the
    # corrections to the mean-field value. We expand to the 2nd order to see
    # the first non-constant term.
    nu_series = nu.series(u, 0, 2)

    # The result of nu_series is `1/2 - c1*u/4 + O(u**2)`.
    # Let's print this in a more readable equation format.
    # The term `-c1*u/4` is the first correction to the mean-field value `1/2`.
    # This correction is of order u^1.
    print("The series expansion for the critical exponent ν(u) around u=0 is:")
    
    # We explicitly print each number in the final equation.
    term_0_num = 1
    term_0_den = 2
    
    term_1_num_coeff = 1
    term_1_den_coeff = 4

    print(f"ν(u) = {term_0_num}/{term_0_den} - (c1/{term_1_den_coeff})*u + O(u^2)")

    print("\nFrom the expansion, we can see the first non-vanishing contribution (correction) to the mean-field value (1/2) is the term proportional to u^1.")

analyze_critical_exponent_nu()
