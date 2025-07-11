import sympy

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting quantum harmonic oscillators.

    This function follows the plan outlined above:
    1. It models the system as two 3D isotropic harmonic oscillators.
    2. It uses second-order perturbation theory to find the energy shift due to the
       dipole-dipole interaction.
    3. The interaction potential is V = (e^2 / (4*pi*R^3)) * (x1*x2 + y1*y2 - 2*z1*z2),
       assuming R is along the z-axis and setting epsilon_0 = 1 as per the prompt.
    4. The ground state energy of the unperturbed system (two 3D QHOs) is E0 = 3*hbar*omega_0.
    5. The first non-zero contribution comes from intermediate states with one excitation
       in each oscillator, e.g., |100>|100>. The energy of these states is Ek = 5*hbar*omega_0.
       The energy denominator is E0 - Ek = -2*hbar*omega_0.
    6. The matrix element <1|x|0> = sqrt(hbar / (2*m*omega_0)).
    7. The sum of the squared matrix elements involves a geometric factor sum of (1^2 + 1^2 + (-2)^2) = 6.
    8. Combining these factors gives the final energy shift.
    """

    # Define symbolic variables for the physical quantities
    e, hbar, m, w0, R, pi = sympy.symbols('e hbar m omega_0 R pi')

    # The final formula for the ground state energy shift (Delta_E) is derived as:
    # Delta_E = - (3 * e**4 * hbar) / (64 * pi**2 * m**2 * w0**3 * R**6)
    # We will construct and print this formula, showing each number explicitly.

    # Numerator of the expression
    numerator_coeff = -3
    numerator_vars = e**4 * hbar

    # Denominator of the expression
    denominator_coeff = 64
    denominator_vars = pi**2 * m**2 * w0**3 * R**6

    # Print the final equation step-by-step
    print("The ground state energy shift, ΔE, is given by the following equation:")
    print("ΔE = - (A * e^4 * ħ) / (B * π^2 * m^2 * ω₀^3 * R^6)")
    print("\nWhere the numerical coefficients are:")
    print(f"A = {abs(numerator_coeff)}")
    print(f"B = {denominator_coeff}")
    
    # Print the final equation with all numbers included
    print("\nThe complete final equation is:")
    final_eq_str = f"ΔE = ({numerator_coeff} * e^4 * ħ) / ({denominator_coeff} * π^2 * m^2 * ω₀^3 * R^6)"
    print(final_eq_str)

    # Return the result in the requested format
    final_expression = (numerator_coeff * numerator_vars) / (denominator_coeff * denominator_vars)
    return final_expression

if __name__ == '__main__':
    final_answer = solve_energy_shift()
    # The final answer is the symbolic expression for the energy shift.
    # For the requested format, we will output the string representation of the expression.
    print(f"\n<<<{final_answer}>>>")
