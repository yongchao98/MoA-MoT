def describe_correlation_transform():
    """
    Explains the space-time Fourier transform of the generalized pair correlation function
    in the context of nuclear criticality.
    """
    name = "Power Spectral Density (PSD)"
    equation_symbol = "S(k, ω)"
    correlation_function = "g(r, τ)"

    print(f"The space-time, double Fourier transform of the generalized pair correlation function ({correlation_function}) is commonly called the '{name}'.")
    print("\n--- Conceptual Explanation ---")
    print(f"The pair correlation function, {correlation_function}, measures the probability of finding a neutron at a certain position and time relative to another neutron.")
    print("Applying a Fourier transform in both space (r -> k, wavevector) and time (τ -> ω, frequency) yields a function in the frequency domain.")
    print(f"This transformed function is denoted as {equation_symbol} and is defined as:")
    print(f"  {equation_symbol} = ∫∫ {correlation_function} * e^(-i(k·r - ωτ)) dr dτ")
    
    print("\n--- Result in a Simplified Model ---")
    print("For many applications, a simplified 'point reactor' model is used. In this case, the spatial dependence is ignored,")
    print("and the PSD's dependence on angular frequency (ω) takes a characteristic shape known as a Lorentzian:")
    
    final_equation = "PSD(ω) ∝ 1 / (ω² + α²)"
    print(f"\nThe final equation is: {final_equation}")

    print("\n--- Components of the Final Equation ---")
    print("  Symbol 'PSD(ω)': The Power Spectral Density as a function of angular frequency ω.")
    print("  Symbol '∝'    : 'is proportional to'.")
    print("  Number '1'      : The numerator of the fraction.")
    print("  Symbol '/'      : The division operator.")
    print("  Symbol 'ω'      : The angular frequency.")
    print("  Operator '^'    : The exponentiation operator (power of).")
    print("  Number '2'      : The exponent for the frequency (ω squared).")
    print("  Operator '+'    : The addition operator.")
    print("  Symbol 'α'      : The prompt neutron decay constant, also known as 'Rossi-alpha'.")
    print("  Operator '^'    : The exponentiation operator (power of).")
    print("  Number '2'      : The exponent for the alpha constant (α squared).")

describe_correlation_transform()
