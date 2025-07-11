import sympy
from sympy import pi, latex

def solve_energy_shift():
    """
    Calculates the ground state energy shift for two interacting quantum harmonic oscillators
    using second-order perturbation theory.
    """
    # Define the symbolic variables for the physical constants and parameters
    e, m, omega0, R, epsilon0, hbar = sympy.symbols('e m omega0 R epsilon0 hbar', real=True, positive=True)

    # 1. The interaction Hamiltonian H' is proportional to a constant C
    # H' = C * (x1*x2 + y1*y2 - 2*z1*z2)
    C = e**2 / (4 * pi * epsilon0 * R**3)

    # 2. The energy denominator in the 2nd order perturbation formula is E0 - En
    # E0 (ground state energy) for two 3D QHOs = (3/2)*hbar*omega0 * 2 = 3*hbar*omega0
    E0 = 3 * hbar * omega0
    # En (energy of excited state) = (5/2)*hbar*omega0 * 2 = 5*hbar*omega0
    En = 5 * hbar * omega0
    energy_denominator = E0 - En

    # 3. The squared matrix element of the position operator is <1|x|0>^2
    # For a 1D QHO, |<1|x|0>|^2 = hbar / (2*m*omega0)
    pos_matrix_element_sq = hbar / (2 * m * omega0)

    # 4. Calculate the contributions from x, y, and z terms in H'
    # The matrix element for the H' term is, e.g., for x: <n|C*x1*x2|0> = C * <1|x1|0> * <1|x2|0>
    # So its square is C^2 * |<1|x1|0>|^2 * |<1|x2|0>|^2 = C^2 * pos_matrix_element_sq^2
    
    # Contribution from the 'x1*x2' term (coefficient is 1)
    term_x = (1 * C)**2 * pos_matrix_element_sq**2 / energy_denominator

    # Contribution from the 'y1*y2' term (coefficient is 1)
    term_y = (1 * C)**2 * pos_matrix_element_sq**2 / energy_denominator

    # Contribution from the '-2*z1*z2' term (coefficient is -2)
    term_z = (-2 * C)**2 * pos_matrix_element_sq**2 / energy_denominator

    # 5. The total energy shift is the sum of these three terms
    Delta_E = term_x + term_y + term_z

    # Simplify the final expression
    Delta_E_simplified = sympy.simplify(Delta_E)

    # 6. Output the final equation and its components
    print("The final equation for the ground state zero energy shift is:")
    # Use sympy.pretty for a cleaner console output
    print(f"\nΔE = {sympy.pretty(Delta_E_simplified, use_unicode=True)}\n")
    
    # Extract coefficients and powers as requested
    num, den = sympy.fraction(Delta_E_simplified)
    
    # Get the main numerical coefficient from the numerator and denominator
    coeff_num = -num.as_coeff_mul()[0] # The expression is negative
    coeff_den = den.as_coeff_mul(pi**2)[0]

    # Use a dictionary to store powers of each symbol for a structured output
    powers = {
        'hbar': sympy.degree(Delta_E_simplified, hbar),
        'e': sympy.degree(Delta_E_simplified, e),
        'epsilon0': -sympy.degree(den, epsilon0),
        'm': -sympy.degree(den, m),
        'omega0': -sympy.degree(den, omega0),
        'pi': -sympy.degree(den, pi),
        'R': -sympy.degree(den, R)
    }

    print("To explicitly output each number in the final equation:")
    print(f"Numerical coefficient in the numerator: {coeff_num}")
    print(f"Numerical coefficient in the denominator: {coeff_den}")
    print(f"Power of Planck's constant (hbar): {powers['hbar']}")
    print(f"Power of elementary charge (e): {powers['e']}")
    print(f"Power of vacuum permittivity (epsilon0): {powers['epsilon0']}")
    print(f"Power of mass (m): {powers['m']}")
    print(f"Power of angular frequency (omega0): {powers['omega0']}")
    print(f"Power of Pi (pi): {powers['pi']}")
    print(f"Power of distance (R): {powers['R']}")
    
    # Return the final expression string as requested for the answer format.
    return Delta_E_simplified

if __name__ == '__main__':
    final_expression = solve_energy_shift()
    # The final answer in the required format
    print(f"\n<<<{final_expression}>>>")
