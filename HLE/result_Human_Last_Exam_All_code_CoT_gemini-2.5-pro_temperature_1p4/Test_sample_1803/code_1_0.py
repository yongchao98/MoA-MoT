import sympy

def calculate_energy_shift():
    """
    Calculates the leading term of the ground state energy shift for two
    interacting quantum harmonic oscillators using second-order perturbation theory.
    """
    # Define symbolic variables for the physical quantities
    e, m, omega0, hbar, R, pi = sympy.symbols('e m omega0 hbar R pi', real=True, positive=True)

    # --- Step 1: Define constants and coefficients from the problem ---

    # The interaction constant from the Coulomb potential, as given e^2/(4*pi*r)
    # We define C_int = e^2 / (4*pi)
    C_int = e**2 / (4 * pi)

    # The coefficient of the dipole-dipole interaction Hamiltonian
    # H_int = C_pert * (x1*x2 + y1*y2 - 2*z1*z2)
    C_pert = C_int / R**3

    # --- Step 2: Define terms for second-order perturbation theory ---

    # Energy denominator: E_0 - E_n = -2*hbar*omega_0
    # where the excited state |n> has both oscillators excited by one quantum.
    energy_denominator = -2 * hbar * omega0

    # Squared matrix element of the position operator between ground and first excited state
    # |<1|x|0>|^2 = hbar / (2*m*omega_0)
    pos_matrix_element_sq = hbar / (2 * m * omega0)

    # --- Step 3: Calculate the contribution from each term in H_int ---

    # The squared matrix element |<n|H_int|0>|^2 for each component
    # For the x-term: |<1x,1x| C_pert * x1*x2 |0,0>|^2 = C_pert^2 * |<1|x|0>|^4
    # Note: <1x,1x|...|0,0> involves two position operators, so we have pos_matrix_element_sq twice.
    mat_el_sq_x = C_pert**2 * pos_matrix_element_sq**2
    mat_el_sq_y = C_pert**2 * pos_matrix_element_sq**2
    mat_el_sq_z = (-2 * C_pert)**2 * pos_matrix_element_sq**2

    # --- Step 4: Sum the contributions to get the total energy shift ---

    # Delta_E = Sum of ( |<n|H_int|0>|^2 / (E0-En) )
    delta_E_x = mat_el_sq_x / energy_denominator
    delta_E_y = mat_el_sq_y / energy_denominator
    delta_E_z = mat_el_sq_z / energy_denominator

    total_energy_shift = sympy.simplify(delta_E_x + delta_E_y + delta_E_z)

    # --- Step 5: Print the final result in a structured format ---

    # Extract the coefficient and powers to display them clearly
    # Expected form: - (3/4) * (e**2/(4*pi))**2 * hbar / (m**2 * omega0**3 * R**6)
    
    # We can programmatically extract parts, but for clarity, we'll format based on our known derivation.
    numerator_coeff = 3
    denominator_coeff = 4
    c_int_power = 2
    c_int_denom = 4
    m_power = 2
    omega_power = 3
    R_power = 6
    
    print("The ground state energy shift is the second-order perturbation correction, Delta E.")
    print("\nThe final calculated expression for the energy shift is:")
    
    # Using an f-string to clearly show each number in the final equation
    final_expression = (
        f"Delta E = - ({numerator_coeff}/{denominator_coeff}) * "
        f"(e**{c_int_power} / ({c_int_denom}*pi))**{c_int_power} * "
        f"hbar / (m**{m_power} * omega0**{omega_power} * R**{R_power})"
    )

    # Re-writing the expression to exactly match the symbolic calculation result for verification
    # result = -3/4 * (e^2/(4*pi*R^3))^2 * hbar / (m^2 * omega0^3)
    # result = -3/4 * e^4/(16*pi^2*R^6) * hbar / (m^2 * omega0^3)
    final_expression_simplified = (
        f"Delta E = -({numerator_coeff * 1}) * hbar * e**{c_int_power*2} / "
        f"({denominator_coeff * c_int_denom**2} * pi**{c_int_power} * m**{m_power} * omega0**{omega_power} * R**{R_power})"
    )

    print("\nIn standard form:")
    final_expression_standard = (
        f"Delta E = - ({numerator_coeff}/{denominator_coeff}) * "
        f"(e**2 / ({c_int_denom}*pi))**{c_int_power} * "
        f"(hbar / (m**{m_power} * omega0**{omega_power} * R**{R_power}))"
    )

    print(final_expression_standard)

    print("\nSymbolic result from sympy:", total_energy_shift)


calculate_energy_shift()
<<<Delta E = - (3/4) * (e**2 / (4*pi))**2 * (hbar / (m**2 * omega0**3 * R**6))>>>