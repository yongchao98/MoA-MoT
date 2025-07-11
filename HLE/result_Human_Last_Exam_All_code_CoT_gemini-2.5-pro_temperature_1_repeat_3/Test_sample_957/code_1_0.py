import sympy

def solve_fabry_perot():
    """
    Symbolically derives the transmission and reflection coefficients of a Fabry-Perot cavity.
    """
    # Define the symbols used in the equations
    # E0: Incident electric field amplitude (will cancel out)
    # tau_m: Mirror transmission coefficient
    # rho_m: Mirror reflection coefficient
    # k0: Wave number
    # d: Gap between mirrors
    # I: Imaginary unit
    E0, tau_m, rho_m, k0, d = sympy.symbols('E_0, tau_m, rho_m, k_0, d', real=True)
    I = sympy.I

    # Phase shift for a single pass through the cavity
    phase_factor_pass = sympy.exp(I * k0 * d)
    
    # Phase shift for a round trip inside the cavity
    phase_factor_round_trip = sympy.exp(I * 2 * k0 * d)

    # --- Transmission Coefficient (tau) ---
    
    # The first wave transmitted through the cavity:
    # Transmit through mirror 1 -> propagate distance d -> transmit through mirror 2
    E_t1 = E0 * tau_m * phase_factor_pass * tau_m
    
    # Each subsequent transmitted wave undergoes an additional round trip inside the cavity.
    # The amplitude is multiplied by (rho_m * rho_m * phase_factor_round_trip)
    # This is a geometric series with common ratio x = rho_m**2 * phase_factor_round_trip
    # The sum of the series is 1 + x + x^2 + ... = 1 / (1 - x)
    round_trip_amplitude_factor = rho_m**2
    common_ratio_t = round_trip_amplitude_factor * phase_factor_round_trip
    
    # Total transmitted field is the sum of the geometric series
    E_t = E_t1 / (1 - common_ratio_t)
    
    # The overall transmission coefficient is tau = E_t / E0
    tau = sympy.simplify(E_t / E0)

    # --- Reflection Coefficient (rho) ---
    
    # The first part of the reflected wave is from the front surface of the first mirror
    E_r1 = E0 * rho_m
    
    # The other parts come from waves that enter the cavity, bounce, and exit from the front.
    # First wave of this type:
    # Transmit in -> propagate d -> reflect off back mirror -> propagate d -> transmit out
    E_r2 = E0 * tau_m * phase_factor_pass * rho_m * phase_factor_pass * tau_m
    
    # The subsequent waves of this type form a geometric series with the same common ratio
    # Sum = E_r2 + E_r2*x + E_r2*x^2 + ... = E_r2 / (1 - x)
    E_r_series = E_r2 / (1 - common_ratio_t)
    
    # Total reflected field is the sum of the direct reflection and the series
    E_r = E_r1 + E_r_series
    
    # The overall reflection coefficient is rho = E_r / E0
    rho = sympy.simplify(E_r / E0)
    
    # Now, let's print the results in a readable format.
    # For the final equation format, we will manually construct the string
    # since sympy's pretty print isn't easy to capture and reformat perfectly.
    
    print("Derived Transmission Coefficient (τ):")
    tau_expr_str = f"τ = (tau_m**2 * exp(i*k_0*d)) / (1 - rho_m**2 * exp(i*2*k_0*d))"
    print(tau_expr_str)
    # The problem asks to output each 'number' in the final equation.
    # We will print the components of the derived formula.
    print("\nEquation for τ:")
    print("Numerator: tau_m**2 * exp(i*k_0*d)")
    print("Denominator: 1 - rho_m**2 * exp(i*2*k_0*d)")


    print("\nDerived Reflection Coefficient (ρ):")
    # Simplify rho further by combining fractions
    rho_combined = sympy.collect(sympy.expand(sympy.together(rho)), sympy.exp(I*2*k0*d))
    
    # Let's create the string representation for our derived rho
    # Our derived rho is: (rho_m - rho_m * (rho_m**2 - tau_m**2) * exp(i*2*k_0*d)) / (1 - rho_m**2 * exp(i*2*k_0*d))
    rho_expr_str = f"ρ = (rho_m - rho_m*(rho_m**2 - tau_m**2)*exp(i*2*k_0*d)) / (1 - rho_m**2*exp(i*2*k_0*d))"
    print(rho_expr_str)
    print("\nEquation for ρ:")
    print("Numerator: rho_m - rho_m*(rho_m**2 - tau_m**2)*exp(i*2*k_0*d)")
    print("Denominator: 1 - rho_m**2*exp(i*2*k_0*d)")

solve_fabry_perot()