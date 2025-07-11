import sympy as sp

def solve_photon_gas_properties():
    """
    Calculates and prints the equilibrium mean energy and entropy of a photon gas
    using symbolic mathematics.
    """
    # Step 1: Define all necessary symbols
    # We define the final quantities E (Energy) and S (Entropy) for display purposes.
    # The variables for calculation are energy (eps), Temperature (T), Volume (V),
    # Boltzmann constant (k_B), speed of light (c), and reduced Planck constant (hbar).
    E, S = sp.symbols('E S')
    eps, T, V, k_B, c, hbar = sp.symbols('epsilon T V k_B c hbar', positive=True)
    
    # The Lagrange multiplier beta is related to temperature
    beta = 1 / (k_B * T)

    # Step 2: Define the key physical formulas
    
    # The density of states g(eps) for photons in a 3D box of volume V.
    # This formula tells us how many available states there are at a given energy.
    # Note: This is often written with h, where h = 2*pi*hbar. We use hbar.
    g_eps = (V / (sp.pi**2 * c**3 * hbar**3)) * eps**2

    # The Bose-Einstein distribution n(eps) for photons (which are bosons with zero chemical potential).
    # This gives the average occupation number of a state with energy eps.
    n_eps = 1 / (sp.exp(eps / (k_B * T)) - 1)

    # Step 3: Calculate the equilibrium mean energy (E)
    # E is the integral from 0 to infinity of energy * density_of_states * occupation_number.
    print("Calculating the mean energy E by integrating eps * g(eps) * n(eps) from 0 to infinity...")
    
    integrand_E = eps * g_eps * n_eps
    
    # Sympy can solve this definite integral.
    # The integral is of the form integral(x**3 / (exp(x) - 1)), which is a known form
    # related to the Riemann zeta function, evaluating to pi**4 / 15.
    mean_E = sp.integrate(integrand_E, (eps, 0, sp.oo))
    
    # Step 4: Calculate the equilibrium entropy (S)
    # For a photon gas (or any gas of massless relativistic particles), there is a
    # simple thermodynamic relation between entropy S, energy E, and temperature T.
    print("Calculating the entropy S using the thermodynamic relation S = (4/3) * E / T...")
    entropy_S = sp.simplify((sp.Rational(4, 3) * mean_E / T))

    # Step 5: Print the final results in a clear format
    # The expressions show all the constants and variables, including the numerical factors.
    print("\n--- Equilibrium Values ---")
    
    # Create formatted strings for the equations
    mean_E_str = f"E = {mean_E}"
    entropy_S_str = f"S = {entropy_S}"
    
    print(f"The equilibrium mean energy is:\n{mean_E_str}\n")
    print(f"The equilibrium entropy is:\n{entropy_S_str}\n")

    # Return the final answer in the specified format for the platform
    final_answer_str = f"E = {str(mean_E)}; S = {str(entropy_S)}"
    return final_answer_str

# Execute the function and capture the final answer string
final_answer = solve_photon_gas_properties()

# The final answer is printed in the required format at the end.
# print(f"<<<{final_answer}>>>")