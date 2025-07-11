import math

def solve_energy_shift():
    """
    This function calculates and prints the symbolic formula for the
    ground state zero-point energy shift between two quantum harmonic oscillators.
    """

    # The interaction is the London dispersion force, with an energy shift
    # proportional to 1/R^6. The general formula derived from second-order
    # perturbation theory is:
    # ΔE = - C / R^6
    # Our goal is to find the expression for the coefficient C.

    # After performing the second-order perturbation calculation, the formula is found to be:
    # ΔE = - (3/4) * (ħ / (m^2 * ω_0^3)) * (e^2 / (4*π*ε_0))^2 * (1 / R^6)

    # Let's define the numerical constants and powers for the final equation.
    
    # Numerical coefficient in the numerator
    coeff_A = 3
    
    # Numerical coefficient in the denominator
    coeff_B = 4
    
    # Power of mass 'm'
    power_m = 2
    
    # Power of frequency 'ω_0'
    power_omega = 3
    
    # Power of the Coulomb term 'e^2 / (4*π*ε_0)'
    power_coulomb = 2
    
    # Power of distance 'R'
    power_R = 6

    print("The leading term for the ground state energy shift (ΔE) is calculated using second-order perturbation theory.")
    print("The final formula is constructed below, with each numerical component explicitly stated.")
    print("-" * 50)
    
    # We construct and print the final equation in a clear, readable format.
    # This fulfills the request to "output each number in the final equation!".
    print("Final Equation for the Energy Shift:")
    print(
        f"ΔE = - ( {coeff_A} * ħ * (e^2 / (4*π*ε_0))^{power_coulomb} ) / "
        f"( {coeff_B} * m^{power_m} * ω_0^{power_omega} * R^{power_R} )"
    )
    print("-" * 50)


# Execute the function to print the result.
solve_energy_shift()