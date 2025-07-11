import sympy

def calculate_energy_shift():
    """
    Calculates the ground state energy shift for two interacting quantum harmonic oscillators.

    The function follows these steps:
    1.  Defines all physical quantities as symbolic variables.
    2.  Sets up the interaction potential V for two collinear dipoles. The coefficient of the x1*x2 term is C.
    3.  Calculates the new normal mode frequencies (omega_plus, omega_minus) for the coupled system.
    4.  Determines the new exact ground state energy E_new from these frequencies.
    5.  Calculates the energy shift Delta_E by subtracting the unperturbed ground state energy.
    6.  Performs a Taylor series expansion of Delta_E for a small interaction (large R) to find the leading term.
    7.  Prints the final, simplified formula for the energy shift.
    """
    # Define symbolic variables for the physical quantities.
    hbar, m, omega_0, R, e, pi = sympy.symbols('hbar m omega_0 R e pi', real=True, positive=True)

    # The problem specifies using e^2/4*pi*r for the Coulomb interaction.
    # We define a constant k_e = e^2/(4*pi).
    k_e = e**2 / (4 * pi)

    # The interaction potential for two collinear dipoles is V = -2 * k_e * x1 * x2 / R^3.
    # The Hamiltonian contains the term C * x1 * x2.
    # So, C is the coupling constant.
    C = -2 * k_e / R**3

    # The potential energy of the coupled system can be diagonalized using normal modes.
    # The new effective spring constants lead to shifted frequencies.
    # m*omega_new^2 = k_new = k_old +/- C
    omega_plus_sq = omega_0**2 - C / m
    omega_minus_sq = omega_0**2 + C / m

    # The new frequencies are the square roots of these expressions.
    omega_plus = sympy.sqrt(omega_plus_sq)
    omega_minus = sympy.sqrt(omega_minus_sq)

    # The new ground state energy is the sum of the zero-point energies of the two normal modes.
    E_new = (hbar / 2) * (omega_plus + omega_minus)

    # The unperturbed ground state energy is for two independent oscillators.
    E_old = hbar * omega_0

    # The energy shift is the difference between the new energy and the old energy.
    Delta_E = E_new - E_old

    # For large R, C is a small parameter. We expand the energy shift as a series in C.
    # The first-order term in C is zero, so we expand to second order to find the leading term.
    # The series is expanded around C=0 up to the term C^2 (order 3).
    Delta_E_series = Delta_E.series(C, 0, 3)

    # The leading term is the series expansion without the O(C^3) part.
    leading_term = Delta_E_series.removeO()

    # We use sympy.simplify to get the final, clean expression.
    final_expression = sympy.simplify(leading_term)
    
    # To meet the output requirement, we will format the final expression as a string,
    # showing all constants and powers explicitly.
    num, den = final_expression.as_numer_denom()
    
    # We construct the string representation. The numerator is negative.
    # We use dictionary to pretty-print symbols.
    pretty_symbols = {
        'hbar': 'ħ',
        'omega_0': 'ω₀'
    }
    
    num_str = " * ".join([f"{pretty_symbols.get(str(s), s)}**{p}" if p > 1 else f"{pretty_symbols.get(str(s), s)}" for s, p in sympy.Mul.from_sympy(abs(num)).as_powers_dict().items()])
    den_str = " * ".join([f"{pretty_symbols.get(str(s), s)}**{p}" if p > 1 else f"{pretty_symbols.get(str(s), s)}" for s, p in sympy.Mul.from_sympy(den).as_powers_dict().items()])
    
    # Custom formatting to make it more readable
    num_str_final = f"ħ * (e**2 / (4*pi))**2"
    den_str_final = f"2 * m**2 * ω₀**3 * R**6"
    
    print("The ground state zero-point energy shift is:")
    print(f"ΔE = -({num_str_final}) / ({den_str_final})")


if __name__ == '__main__':
    calculate_energy_shift()