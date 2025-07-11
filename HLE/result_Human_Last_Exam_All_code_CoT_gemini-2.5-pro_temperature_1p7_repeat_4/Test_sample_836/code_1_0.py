import sympy as sp

def define_dynamic_structure_factor():
    """
    Symbolically defines and prints the equation for the Dynamic Structure Factor.

    The Dynamic Structure Factor S(k, omega) is the space-time, double Fourier
    transform of the generalized pair correlation function G(r, t). This function
    is crucial in analyzing particle fluctuations, such as in neutron noise
    analysis for nuclear criticality.
    """

    # Define the symbols
    # G(r, t): The generalized pair correlation function (Van Hove function)
    # S(k, omega): The Dynamic Structure Factor
    r, t, k, omega = sp.symbols('r t k omega', real=True)
    G = sp.Function('G')(r, t)
    S = sp.Function('S')(k, omega)
    I = sp.I # Imaginary unit

    # Define the Fourier transform kernel
    # The sign convention exp(i*(omega*t - k*r)) is common in this field
    fourier_kernel = sp.exp(I * (omega * t - k * r))

    # Define the integrand
    integrand = G * fourier_kernel

    # Create the integral expression for the double Fourier transform
    # The integration is over all space (r) and all time (t)
    transform_integral = sp.Integral(integrand, (r, -sp.oo, sp.oo), (t, -sp.oo, sp.oo))

    # Create the full equation
    equation = sp.Eq(S, transform_integral)

    # Print the explanation and the symbolic equation
    print("The space-time, double Fourier transform of the generalized pair correlation function G(r, t)")
    print("is known as the Dynamic Structure Factor, S(k, omega).")
    print("\nIts definition is given by the following equation (shown for one spatial dimension 'r'):\n")

    # The sp.pretty_print function provides a nice text-based representation of the equation.
    sp.pretty_print(equation)

    print("\nWhere:")
    print("  S(k, \u03C9) is the Dynamic Structure Factor")
    print("  G(r, t) is the generalized pair correlation function")
    print("  k is the wave vector (spatial frequency)")
    print("  \u03C9 (omega) is the angular frequency")
    print("  r is the spatial separation")
    print("  t is the time separation")
    print("  i is the imaginary unit")

if __name__ == '__main__':
    define_dynamic_structure_factor()
