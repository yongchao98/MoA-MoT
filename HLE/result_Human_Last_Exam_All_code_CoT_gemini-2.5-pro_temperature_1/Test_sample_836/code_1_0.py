import sympy

def describe_fourier_transform():
    """
    Symbolically describes the space-time Fourier transform of the
    generalized pair correlation function and identifies its common name in
    nuclear criticality.
    """

    # Define symbolic variables
    # r represents the spatial separation vector between two detection events.
    # t represents the time lag between two detection events.
    # k represents the wave vector, the Fourier conjugate to space r.
    # omega represents the angular frequency, the Fourier conjugate to time t.
    r, t, k, omega = sympy.symbols('r t k omega')

    # Define the functions
    # C(r, t) is the generalized pair correlation function in the space-time domain.
    C_rt = sympy.Function('C')(r, t)

    # S(k, omega) is the result of the double Fourier transform.
    S_kw = sympy.Function('S')(k, omega)

    # The Fourier Transform operation is a conceptual representation here.
    # The equation is: S(k, ω) = ∬ C(r, t) * e^(-i(k·r - ωt)) dr dt
    # We will represent this symbolically.
    
    print("In nuclear criticality, the relationship between neutron detection events is described by the generalized pair correlation function.")
    print(f"This function is typically written as: {C_rt}")
    print("\nApplying a double Fourier transform converts this function from the space-time domain (r, t) to the wavevector-frequency domain (k, \u03C9).") # \u03C9 is unicode for ω
    
    # Create the symbolic equation
    # We use a custom 'FT' function to represent the Fourier Transform operator for clarity.
    FT_op = sympy.Function('FT')
    symbolic_equation = sympy.Eq(S_kw, FT_op(C_rt))

    # Print the final equation with each symbol
    # The 'args' attribute of the symbolic equation lets us access each part.
    # Left side: S(k, ω)
    # Right side: FT(C(r, t))
    print("\nThe transformation is represented by the equation:")
    print(f"{symbolic_equation.lhs.func}({symbolic_equation.lhs.args[0]}, {symbolic_equation.lhs.args[1]}) = {symbolic_equation.rhs.func}[{symbolic_equation.rhs.args[0].func}({symbolic_equation.rhs.args[0].args[0]}, {symbolic_equation.rhs.args[0].args[1]})]")

    print("\nThis resulting function, S(k, \u03C9), is commonly called the Power Spectral Density (PSD).")
    print("It describes the 'power' of the neutron population fluctuations at a specific spatial wavevector (k) and temporal frequency (\u03C9).")


if __name__ == '__main__':
    describe_fourier_transform()
