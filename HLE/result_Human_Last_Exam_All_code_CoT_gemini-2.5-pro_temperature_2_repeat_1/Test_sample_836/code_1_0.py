import sympy

def display_fourier_transform_identity():
    """
    Explains and displays the space-time, double Fourier transform of the
    generalized pair correlation function and provides its common name in the
    nuclear criticality community.
    """
    
    # Define mathematical symbols for clarity
    S_kw = sympy.Symbol("S(k, w)")
    C_rtau = sympy.Symbol("C(r, tau)")
    i = sympy.I
    k, r = sympy.symbols('k, r', real=True, vector=True)
    w, tau = sympy.symbols('omega, tau', real=True)
    
    # Explain the concepts
    print("The 'generalized pair correlation function' in nuclear systems is the covariance of the neutron flux.")
    print("Let's denote it as C(r, tau), where 'r' is the spatial separation and 'tau' is the time lag.\n")
    
    print("The space-time, double Fourier transform of this function is a new function, S(k, w).")
    print("Here, 'k' is the wavevector (the Fourier conjugate of space 'r') and 'w' (omega) is the angular frequency (the Fourier conjugate of time 'tau').\n")
    
    # Display the defining equation
    print("The definitional equation for this transform is:\n")
    
    # Manually print the equation to fulfill the user request
    # to output each 'number' (term) in the equation.
    print("S(k, w)  =  Integral over all tau [ Integral over all r [ C(r, tau) * exp(-i * (k . r - w * tau)) ] d_r ] d_tau")
    print("\nWhere:")
    print(" S(k, w)  => The resulting function in the frequency domain.")
    print(" Integral => A double integral, one over all space (r) and one over all time (tau).")
    print(" C(r, tau)=> The generalized pair correlation function (covariance of neutron flux).")
    print(" exp(...) => The complex exponential which forms the kernel of the Fourier transform.")
    print(" i        => The imaginary unit.")
    print(" k . r    => The dot product of the wavevector and the position vector.")
    print(" w * tau  => The product of angular frequency and time lag.\n")

    # State the final answer
    print("-" * 60)
    print("This function, S(k, w), is most commonly called the:")
    print("\033[1mNoise Spectral Density\033[0m")
    print("Sometimes, it is also referred to more specifically as the Power Spectral Density (PSD)")
    print("or Cross-Power Spectral Density (CPSD) of the neutron noise.")
    print("-" * 60)


if __name__ == '__main__':
    display_fourier_transform_identity()
    print("\n<<<Noise Spectral Density>>>")
