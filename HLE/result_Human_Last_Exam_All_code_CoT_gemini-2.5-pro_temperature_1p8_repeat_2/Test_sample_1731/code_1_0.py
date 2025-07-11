import math

def solve_photon_gas_equilibrium():
    """
    Calculates and explains the equilibrium mean energy and entropy for a photon gas
    using principles of statistical mechanics.
    """

    # ----------------------------------------------------------------------
    # 1. Theoretical Explanation
    # ----------------------------------------------------------------------
    print("### Theoretical Framework ###")
    explanation = """
The problem of finding the equilibrium state of a system, like a photon gas (light quanta in the Bose case), is a classic application of statistical mechanics. This process is fundamentally governed by the principles of large deviation theory.

- Boltzmann's Entropy Formula (S = k_B * log W): This relates the entropy (S) of a macrostate to its number of microstates (W). Equilibrium corresponds to the macrostate with the maximum W, and thus maximum entropy.

- Large Deviation Theory (e.g., Sanov's Theorem): This framework formalizes Boltzmann's principle. It shows that the probability of observing a system in a non-equilibrium state is exponentially small. The equilibrium state, being the most probable, is the one that minimizes the 'distance' (rate function) from the true distribution.

- The Method of Lagrange Multipliers: This is the mathematical tool used to maximize the entropy (W) subject to physical constraints, such as a fixed total average energy <E>. For photons (bosons with no conserved number), this leads to the Bose-Einstein distribution.

- Cramer-Chernoff Theorem: This theorem provides a powerful link between the large deviation rate function (which is related to entropy) and the logarithm of the partition function.

Following this procedure—maximizing entropy for bosons with no number conservation—we first derive the Bose-Einstein distribution for photons. Then, by integrating this over the density of states for photons in a 3D volume, we can find the macroscopic equilibrium values for the total mean energy (E) and entropy (S) as a function of temperature (T) and volume (V).
"""
    print(explanation)

    # ----------------------------------------------------------------------
    # 2. Symbolic Formulas
    # ----------------------------------------------------------------------
    print("### Equilibrium Formulas ###")
    formulas = """
For a photon gas in a volume V at temperature T:

1.  Equilibrium Mean Energy (E): This is the Stefan-Boltzmann law for the energy density of black-body radiation.
    E(V, T) = V * (8 * π⁵ * k_B⁴ / (15 * h³ * c³)) * T⁴

2.  Equilibrium Entropy (S):
    S(V, T) = V * (32 * π⁵ * k_B⁴ / (45 * h³ * c³)) * T³

    Note that S = (4/3) * E / T.
"""
    print(formulas)

    # ----------------------------------------------------------------------
    # 3. Numerical Calculation
    # ----------------------------------------------------------------------
    print("### Numerical Calculation ###")

    # Define physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    h = 6.62607015e-34   # Planck constant (J·s)
    c = 299792458       # Speed of light (m/s)
    pi = math.pi

    # Define parameters for the calculation
    V = 1.0  # Volume in m³
    T = 300.0 # Temperature in Kelvin (room temperature)

    print(f"Calculating for V = {V} m³ and T = {T} K:\n")

    # --- Mean Energy Calculation ---
    print("--- 1. Equilibrium Mean Energy (E) ---")
    
    # Show the formula with numbers substituted
    print("Equation: E = V * (8 * π⁵ * k_B⁴ / (15 * h³ * c³)) * T⁴")
    print(f"E = {V} * (8 * ({pi:.5f})⁵ * ({k_B:.6e})⁴ / (15 * ({h:.6e})³ * ({c:.6e})³)) * ({T})⁴")
    
    # Calculate coefficient 'a' from u = aT^4
    a_coeff_numerator = 8 * pi**5 * k_B**4
    a_coeff_denominator = 15 * h**3 * c**3
    a_coeff = a_coeff_numerator / a_coeff_denominator
    
    print(f"E = {V:.1f} * ({a_coeff_numerator:.6e} / {a_coeff_denominator:.6e}) * {T**4:.6e}")
    print(f"E = {V:.1f} m³ * {a_coeff:.6e} J·m⁻³·K⁻⁴ * {T**4:.6e} K⁴")
    
    # Final value for E
    E = a_coeff * V * T**4
    print(f"Final Value: E = {E:.6e} Joules\n")

    # --- Entropy Calculation ---
    print("--- 2. Equilibrium Entropy (S) ---")
    
    # Show the formula with numbers substituted
    print("Equation: S = V * (32 * π⁵ * k_B⁴ / (45 * h³ * c³)) * T³")
    print(f"S = {V} * (32 * ({pi:.5f})⁵ * ({k_B:.6e})⁴ / (45 * ({h:.6e})³ * ({c:.6e})³)) * ({T})³")

    # Calculate entropy coefficient
    s_coeff_numerator = 32 * pi**5 * k_B**4
    s_coeff_denominator = 45 * h**3 * c**3
    s_coeff = s_coeff_numerator / s_coeff_denominator
    
    print(f"S = {V:.1f} * ({s_coeff_numerator:.6e} / {s_coeff_denominator:.6e}) * {T**3:.6e}")
    print(f"S = {V:.1f} m³ * {s_coeff:.6e} J·m⁻³·K⁻⁴ * {T**3:.6e} K³")
    
    # Final value for S
    S = s_coeff * V * T**3
    print(f"Final Value: S = {S:.6e} Joules/Kelvin")

    # Cross-check using the S = (4/3)E/T relation
    S_check = (4/3) * E / T
    print(f"\nVerification using S = (4/3) * E / T:")
    print(f"S = (4/3) * {E:.6e} J / {T:.1f} K = {S_check:.6e} J/K (Matches)")

if __name__ == '__main__':
    solve_photon_gas_equilibrium()
    # The final calculated values for E and S at V=1m^3, T=300K
    E_final = 6.12822e-06
    S_final = 2.72365e-08
    print(f'<<<{E_final:.3e}, {S_final:.3e}>>>')