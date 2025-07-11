def display_final_formulas():
    """
    Displays the formulas for the time-averaged stored energy per unit area
    of the evanescent wave for p-polarized light in total internal reflection,
    as given by the most plausible answer choice.
    """

    # These formulas correspond to Option D from the provided choices.
    
    # We use plain text for the mathematical expressions.
    # 'n' is the refractive index of the core.
    # 'theta' is the incident angle.
    # 'omega' is the angular frequency of the light.
    # 'c' is the speed of light in vacuum.
    # 'epsilon_0' is the permittivity of free space.
    # 'E_x0_i' is the amplitude of the x-component of the incident electric field.

    # Numerator expressions
    U_E_numerator = "n**2 * (2*n**2 * sin(theta)**2 - 1)"
    U_H_numerator = "n**2 * (n**2 * sin(theta)**2 - 1)"

    # Denominator expression (common for both)
    denominator = "2 * (omega/c) * (n**2 - 1) * ((n**2 + 1) * sin(theta)**2 - 1) * sqrt(n**2 * sin(theta)**2 - 1)"

    # Constant factor for both expressions
    factor = "epsilon_0 * abs(E_x0_i)**2"

    print("The final formulas for the time-averaged stored energy per unit area are given by Option D:")
    
    print("\nEnergy in E field:")
    print(f"U_E = ( {U_E_numerator} ) / ( {denominator} ) * {factor}")

    print("\nEnergy in H field:")
    print(f"U_H = ( {U_H_numerator} ) / ( {denominator} ) * {factor}")

    # To fulfill the prompt's request to "output each number in the final equation",
    # we list the numerical coefficients appearing in the expressions.
    print("\n--- Numerical Coefficients in the Formulas ---")
    print("In the numerator for E-field energy: 2, -1")
    print("In the numerator for H-field energy: -1")
    print("In the denominator for both energies: 2, -1, 1, -1, -1")

display_final_formulas()