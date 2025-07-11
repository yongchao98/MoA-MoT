import sympy as sp

def solve_photon_gas_equilibrium():
    """
    Calculates and prints the equilibrium values of mean energy and entropy
    for a Bose gas of light quanta (photons) using symbolic mathematics.
    
    The method is based on statistical mechanics, whose principle of maximizing
    entropy is justified by large deviation theorems (Boltzmann-Sanov).
    """

    # --- Step 1: Define all symbolic variables ---
    # We treat physical constants and system parameters as symbols to get an exact formula.
    E, S = sp.symbols('E S', real=True)  # Target variables: Mean Energy, Entropy
    e, T, V = sp.symbols('e T V', positive=True, real=True) # State variables: energy level, Temperature, Volume
    k_B = sp.Symbol('k_B', positive=True, real=True, constant=True) # Boltzmann constant
    h_bar = sp.Symbol('hbar', positive=True, real=True, constant=True) # Reduced Planck constant
    c = sp.Symbol('c', positive=True, real=True, constant=True) # Speed of light

    print("Deriving equilibrium values for a photon gas...")
    print("The calculation uses the Bose-Einstein distribution for photons and the density of states in 3D.")
    
    # --- Step 2: Define the distribution and density of states ---
    
    # The Bose-Einstein distribution for the average occupation number <n(e)> of a state with energy 'e'.
    # For photons, the chemical potential is zero as their number is not conserved.
    # beta = 1 / (k_B * T)
    n_e = 1 / (sp.exp(e / (k_B * T)) - 1)

    # The density of states g(e) for photons (spin-1, two polarization states) in a 3D volume V.
    g_e = (V / (sp.pi**2 * h_bar**3 * c**3)) * e**2

    # --- Step 3: Calculate the Equilibrium Mean Energy (E) ---
    # The total energy E is the integral of (energy per state) * (occupation number) * (density of states)
    # over all possible energies.
    # E = Integral[e * n(e) * g(e) de] from 0 to infinity.
    
    integrand_E = e * n_e * g_e
    
    # To solve this integral, we perform a change of variables: x = e / (k_B * T)
    # This transforms the integral into a standard form known as the Bose-Einstein integral.
    x = sp.symbols('x', positive=True)
    
    # The integral part is of the form Integral[x^3 / (exp(x) - 1) dx] from 0 to oo.
    # This is related to the Riemann zeta function, zeta(4), and the Gamma function, Gamma(4).
    # The value of the integral is Gamma(4) * zeta(4) = 6 * (pi^4 / 90) = pi^4 / 15.
    integral_value = sp.integrate(x**3 / (sp.exp(x) - 1), (x, 0, sp.oo))

    # The full expression for energy is the coefficient part multiplied by this integral value.
    # Coefficient = (V / (pi^2 * h_bar^3 * c^3)) * (k_B * T)^4
    energy_coefficient = (V / (sp.pi**2 * h_bar**3 * c**3)) * (k_B * T)**4
    
    mean_energy_expr = energy_coefficient * integral_value
    
    # --- Step 4: Calculate the Equilibrium Entropy (S) ---
    # For a photon gas, there is a simple thermodynamic relation between energy and entropy:
    # Pressure P = E / (3V), and Helmholtz Free Energy F = E - TS = -PV.
    # This leads to E - TS = -E/3, which simplifies to TS = (4/3)E or S = (4/3) * E / T.
    
    entropy_expr = sp.simplify((4 * mean_energy_expr) / (3 * T))

    # --- Step 5: Print the final results ---
    
    print("\n--- Equilibrium Mean Energy (E) ---")
    print("This is the Stefan-Boltzmann law for the total energy of black-body radiation in a cavity of volume V.")
    # The sp.pretty() function provides a nicely formatted text output for the equation.
    final_E_eq = sp.Eq(E, mean_energy_expr)
    print("The equation for the mean energy E is:")
    print(sp.pretty(final_E_eq, use_unicode=True))


    print("\n--- Equilibrium Entropy (S) ---")
    print("The entropy is derived from the energy using the thermodynamic relation S = (4/3)E/T.")
    final_S_eq = sp.Eq(S, entropy_expr)
    print("The equation for the entropy S is:")
    print(sp.pretty(final_S_eq, use_unicode=True))


# Execute the function
solve_photon_gas_equilibrium()
final_energy_expression = "(pi**2 * k_B**4 * T**4 * V) / (15 * hbar**3 * c**3)"
final_entropy_expression = "(4 * pi**2 * k_B**4 * T**3 * V) / (45 * hbar**3 * c**3)"
<<<f"Energy: {final_energy_expression}, Entropy: {final_entropy_expression}">>>