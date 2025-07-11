def solve_energy_shift():
    """
    Calculates and prints the formula for the ground state zero-point energy shift
    of two interacting quantum harmonic oscillators.

    The derivation uses second-order perturbation theory for two 3D isotropic
    quantum harmonic oscillators interacting via a dipole-dipole potential.

    The final formula for the energy shift (Delta_E) is printed, followed by
    an explicit list of all the numbers (coefficients and exponents) in that formula.
    """

    # The derived formula for the ground state energy shift
    formula = "ΔE = - (3 * ħ * e⁴) / (4 * m² * ω₀³ * (4 * π * ε₀)² * R⁶)"

    print("The leading term for the ground state energy shift is given by the following formula:")
    print(formula)
    print("\nHere are the numbers that appear in the final equation:")

    # Explicitly list each number as requested
    print("Numerical coefficient in the numerator: 3")
    print("Numerical coefficient in the denominator: 4")
    print("Exponent for mass 'm': 2")
    print("Exponent for frequency 'ω₀': 3")
    print("Numerical coefficient for pi 'π' inside the parenthesis: 4")
    print("Exponent for the parenthesis '(4 * π * ε₀)': 2")
    print("Exponent for distance 'R': 6")

solve_energy_shift()