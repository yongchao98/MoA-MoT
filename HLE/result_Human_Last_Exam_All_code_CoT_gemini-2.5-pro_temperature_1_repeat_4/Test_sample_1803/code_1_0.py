def calculate_energy_shift():
    """
    This function calculates and prints the formula for the leading term of the ground state 
    energy shift for two interacting quantum harmonic oscillators.

    The derivation involves second-order perturbation theory, yielding the van der Waals 
    (London dispersion) interaction energy, which is proportional to 1/R^6.
    """

    # Symbolic representations of the physical quantities
    delta_E = "ΔE"
    h_bar = "ħ"
    e = "e"
    m = "m"
    omega_0 = "ω₀"
    R = "R"
    pi = "π"

    # Numerical coefficients and powers derived from the calculation
    # The final formula is of the form: - (coeff_num / coeff_den) * (h_bar * e^4) / (pi^2 * m^2 * omega_0^3 * R^6)
    
    # Numerator coefficient from perturbation theory calculation
    numerator_coefficient = 3
    
    # Denominator coefficient from perturbation theory calculation ( (4*pi)^2 * 4 = 64 )
    denominator_coefficient = 64
    
    # Powers of the physical quantities in the final expression
    power_of_e = 4
    power_of_pi = 2
    power_of_m = 2
    power_of_omega_0 = 3
    power_of_R = 6

    # Construct the final equation string for display
    equation = (
        f"{delta_E} = - ( {numerator_coefficient} * {h_bar} * {e}^{power_of_e} ) / "
        f"( {denominator_coefficient} * {pi}^{power_of_pi} * {m}^{power_of_m} * {omega_0}^{power_of_omega_0} * {R}^{power_of_R} )"
    )

    print("The leading term in R for the ground state zero-point energy shift is:")
    print(equation)
    
    # Print each numerical component of the equation as requested
    print("\nThe numerical components of this equation are:")
    print(f"The coefficient in the numerator is: {numerator_coefficient}")
    print(f"The coefficient in the denominator is: {denominator_coefficient}")
    print(f"The power of the elementary charge 'e' is: {power_of_e}")
    print(f"The power of pi 'π' is: {power_of_pi}")
    print(f"The power of the mass 'm' is: {power_of_m}")
    print(f"The power of the frequency 'ω₀' is: {power_of_omega_0}")
    print(f"The power of the distance 'R' is: {power_of_R}")

# Execute the function to display the results
calculate_energy_shift()