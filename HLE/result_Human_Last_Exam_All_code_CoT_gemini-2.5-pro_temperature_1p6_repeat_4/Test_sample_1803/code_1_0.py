import sympy

def calculate_energy_shift():
    """
    Symbolically calculates the ground state energy shift for two interacting
    quantum harmonic oscillators using second-order perturbation theory.
    """
    # 1. Define symbols for the physical quantities.
    # We assume all are positive real numbers.
    e, m, omega0, R, hbar, epsilon0 = sympy.symbols('e m omega0 R hbar epsilon0', real=True, positive=True)

    # 2. Define the interaction Hamiltonian H'.
    # This is the dipole-dipole interaction potential for oscillators oriented
    # parallel to each other and perpendicular to the separation vector.
    # The term 'x₁x₂' represents the product of the displacement coordinates.
    # H' = C * x₁ * x₂
    C = e**2 / (4 * sympy.pi * epsilon0 * R**3)

    # 3. Calculate the necessary terms for second-order perturbation theory.
    # The formula is: ΔE = |<1,1|H'|0,0>|^2 / (E_ground - E_excited)

    # 4. The matrix element <1|x|0> for a single QHO.
    x01 = sympy.sqrt(hbar / (2 * m * omega0))

    # 5. The matrix element of the perturbation H' between the ground state |0,0>
    # and the relevant excited state |1,1>.
    # <1,1|H'|0,0> = C * <1|x₁|0> * <1|x₂|0>
    H_prime_element = C * x01**2

    # 6. The energy difference between the ground state and the excited state.
    # E_ground = E₀₀ = ħω₀
    # E_excited = E₁₁ = 3ħω₀
    energy_denominator = (hbar * omega0) - (3 * hbar * omega0)

    # 7. Calculate the second-order energy shift.
    delta_E = (H_prime_element**2) / energy_denominator

    # 8. Simplify the final expression to get the result.
    final_expression = sympy.simplify(delta_E)

    # 9. Output the final equation and its constituent numbers as requested.
    print("The calculated ground state energy shift ΔE due to the interaction is:")
    # sympy.pprint is used for a more readable mathematical output
    sympy.pprint(final_expression, use_unicode=True)

    print("\n" + "="*60)
    print("Breakdown of the numbers and powers in the final equation:")
    print("="*60)

    # Separate the numerical coefficient from the symbolic part
    coeff, terms = final_expression.as_coeff_Mul()
    print(f"Numerical Coefficient: {coeff}")

    # Get the powers of each symbol in the expression
    powers_dict = terms.as_powers_dict()
    for base, exponent in powers_dict.items():
        print(f"Power of '{base}': {exponent}")
    
    # Return a string representation for the final answer block.
    return str(final_expression)


if __name__ == '__main__':
    final_answer_formula = calculate_energy_shift()
    # The final answer will be captured in the required format after this block.
    # Example format: print(f'<<<{final_answer_formula}>>>')

<<<_placeholder_for_final_answer_>>>