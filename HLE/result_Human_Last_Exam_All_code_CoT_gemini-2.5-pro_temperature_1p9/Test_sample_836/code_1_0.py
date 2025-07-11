import textwrap

def print_psd_explanation_and_equation():
    """
    This function explains and prints the equation for the Power Spectral Density (PSD)
    in the context of nuclear reactor noise analysis.
    """

    # --- Step 1 & 2: Explanation ---
    intro = """
    In the nuclear criticality community, the space-time, double Fourier 
    transform of the generalized pair correlation function is called the 
    Power Spectral Density (PSD), often denoted as S(k, ω).

    When simplifying to a point reactor model (ignoring spatial dependence), 
    this becomes a function of angular frequency (ω) only, S(ω). This 
    function describes the frequency content of the neutron population's 
    statistical fluctuations in a nuclear system. It is fundamental to 
    "neutron noise analysis," which is used to measure subcritical reactivity 
    and other important reactor parameters.
    """
    print(textwrap.dedent(intro).strip())
    print("-" * 50)

    # --- Step 3 & 4: Representative Equation (Point Kinetics) ---
    print("The famous equation for the auto-Power Spectral Density in a simple\n"
          "point kinetics model has a Lorentzian shape and is given by:\n")
    
    # Define the components of the equation as strings
    symbol_psd = "S(ω)"
    equals = "="
    numerator_source_term = "A"
    denominator_alpha_squared = "α²"
    plus = "+"
    denominator_omega_squared = "ω²"
    
    # The constraint asks to "output each number in the final equation".
    # Since there are no numbers, we will output each symbol.
    print(f"The final equation is constructed from the following parts:")
    print(f"  {symbol_psd} {equals} {numerator_source_term} / ({denominator_alpha_squared} {plus} {denominator_omega_squared})\n")
    
    print("Where:")
    print("  S(ω) = The Power Spectral Density as a function of angular frequency.")
    print("    A  = A source term, proportional to the fission rate, detector efficiency,")
    print("         and the Diven factor (variance in neutrons per fission).")
    print("    α  = The prompt neutron decay constant, also known as 'Rossi-alpha'.")
    print("    ω  = The angular frequency (2 * pi * f).")

# Execute the function to print the explanation and the equation.
print_psd_explanation_and_equation()
