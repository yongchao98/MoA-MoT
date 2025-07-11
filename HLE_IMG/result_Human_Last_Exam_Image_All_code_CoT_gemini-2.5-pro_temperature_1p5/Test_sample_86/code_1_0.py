def solve():
    """
    This function prints the energy spectrum for a harmonic oscillator with a quartic
    perturbation, as predicted by the tadpole self-energy diagram.
    """
    
    # The derived renormalized frequency squared
    omega_prime_squared_str = "omega_0**2 + (u * hbar**2) / (4 * m**2 * omega_0)"

    # The final energy spectrum formula
    equation_str = f"E_n = hbar * sqrt({omega_prime_squared_str}) * (n + 1/2)"

    print("The tadpole self-energy diagram resums a class of interactions,")
    print("leading to a renormalized harmonic oscillator frequency.")
    print("The predicted energy spectrum E_n is given by the formula:")
    print("")
    print(equation_str)
    print("")
    
    print("In this equation, the numerical constants and exponents are:")
    print("In the term (n + 1/2): numerator is 1, denominator is 2.")
    print("In the square root, which is an exponent of 1/2: numerator is 1, denominator is 2.")
    print("In the correction term (u * hbar**2) / (4 * m**2 * omega_0):")
    print("  - The exponent of hbar is: 2")
    print("  - The numerical coefficient in the denominator is: 4")
    print("  - The exponent of m is: 2")

solve()
