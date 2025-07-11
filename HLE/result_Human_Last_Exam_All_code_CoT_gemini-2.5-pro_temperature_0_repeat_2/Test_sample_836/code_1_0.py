import sympy

def explain_fourier_transform_of_correlation_function():
    """
    This script defines and explains the space-time Fourier transform
    of the generalized pair correlation function in nuclear physics.
    """

    # --- 1. Define the mathematical symbols using SymPy ---
    # We use 1D space for simplicity, but the concept extends to 3D vectors.
    r, tau = sympy.symbols('r, tau', real=True)
    k, omega = sympy.symbols('k, omega', real=True)
    E1, E2 = sympy.symbols('E1, E2', positive=True)
    i = sympy.I

    # Define the functions
    # g: the relative pair correlation function
    g = sympy.Function('g')
    # P: the Power Spectral Density
    P = sympy.Function('P')

    # --- 2. Print the explanation and the result ---
    print("The space-time, double Fourier transform of the generalized pair correlation function")
    print("is commonly called the Power Spectral Density (PSD) of the neutron noise.")
    print("This is a central quantity in nuclear reactor noise analysis and is used to")
    print("infer properties of a nuclear system, such as its reactivity.")
    print("\nUnder the common assumptions of a homogeneous and stationary medium, the correlation")
    print("function g(r1, t1; r2, t2) depends only on the relative separation r = r1 - r2 and")
    print("time lag tau = t1 - t2. Its Fourier transform is the PSD, P(k, omega). This relationship")
    print("is an instance of the Wiener-Khinchin Theorem.")
    print("\n--- The Defining Equation ---")

    # To satisfy the "output each number in the final equation" requirement,
    # we will print the equation in a structured way, explaining each component.

    # Left Hand Side
    lhs = "P(k, E1, omega, E2)"
    # Right Hand Side
    rhs = "Integral(g(r, E1, E2, tau) * exp(-i*k*r) * exp(i*omega*tau), dr, dtau)"

    print(f"{lhs} = {rhs}")
    print("\n--- Components of the Equation ---")
    print(f"  P(k, E1, omega, E2): The Power Spectral Density function.")
    print(f"  g(r, E1, E2, tau)  : The relative pair correlation function.")
    print(f"  Integral(...)      : Integration over all space 'r' and time lag 'tau'.")
    print(f"  exp(...)           : The exponential function, e^x.")
    print(f"  k                  : The wavevector, the Fourier variable conjugate to space 'r'.")
    print(f"  omega              : The angular frequency, the Fourier variable conjugate to time 'tau'.")
    print(f"  E1, E2             : The energies of the two correlated particles (neutrons).")
    print(f"  i                  : The imaginary unit, sqrt(-1).")
    print(f"  r                  : The relative position vector.")
    print(f"  tau                : The time lag.")

if __name__ == '__main__':
    explain_fourier_transform_of_correlation_function()