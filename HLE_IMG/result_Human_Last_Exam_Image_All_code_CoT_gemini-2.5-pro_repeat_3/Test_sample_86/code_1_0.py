import sympy

def solve_energy_spectrum():
    """
    Calculates and prints the energy spectrum for a harmonic oscillator with a quartic perturbation,
    as predicted by the first-order self-energy (tadpole) diagram.
    """

    # Define symbolic variables for clarity, although we will build the final expression as a string.
    hbar = sympy.Symbol('hbar')
    omega0 = sympy.Symbol('omega_0')
    u = sympy.Symbol('u')
    m = sympy.Symbol('m')
    n = sympy.Symbol('n')

    # The self-energy diagram leads to a correction in the squared frequency.
    # The new squared frequency is omega_new^2 = omega_0^2 + delta_omega^2.

    # Step 1: The shift is delta_omega^2 = (u / (2*m)) * <x^2>_0
    # Step 2: The ground state expectation value is <x^2>_0 = hbar / (2 * m * omega_0)
    # Step 3: Combining these, delta_omega^2 = (u / (2*m)) * (hbar / (2*m*omega_0)) = u*hbar / (4*m^2*omega_0)

    # The problem requires printing each number in the final equation.
    # We will construct the final formula as a string.

    # Coefficients
    c1 = 1
    c2 = 2
    c4 = 4

    # Symbolic parts of the equation
    s_En = "E_n"
    s_hbar = "hbar"
    s_n = "n"
    s_omega0_sq = "omega_0**2"
    s_u = "u"
    s_m_sq = "m**2"
    s_omega0 = "omega_0"

    # Building the final equation string
    term_n = f"({s_n} + {c1}/{c2})"
    numerator_corr = f"{c1} * {s_u} * {s_hbar}"
    denominator_corr = f"({c4} * {s_m_sq} * {s_omega0})"
    correction_term = f"({numerator_corr} / {denominator_corr})"
    sqrt_term = f"sqrt({s_omega0_sq} + {correction_term})"
    
    final_equation = f"{s_En} = {c1} * {s_hbar} * {term_n} * {sqrt_term}"

    print("The energy spectrum predicted by the self-energy diagram is given by the formula for a harmonic oscillator with a renormalized frequency.")
    print("The final equation is:")
    print(final_equation)

solve_energy_spectrum()
<<<E_n = 1 * hbar * (n + 1/2) * sqrt(omega_0**2 + ((1 * u * hbar) / (4 * m**2 * omega_0)))>>>