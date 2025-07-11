import math

def derive_bose_equilibrium():
    """
    Derives and prints the equilibrium values for mean energy (U) and entropy (S)
    for a photon gas using principles of statistical mechanics related to
    large deviation theory.
    """

    print("### Derivation of Equilibrium Properties for a Photon Gas (Bose Case) ###")
    print("-" * 70)

    # --- Step 1: Statistical Entropy (Boltzmann-Sanov Principle) ---
    print("1. Statistical Entropy (Connection to Boltzmann-Sanov Theorem)")
    print("The equilibrium state of a system is the macrostate with the maximum number of microstates (W).")
    print("This is the core idea of Boltzmann's entropy formula, S = k_B * ln(W), and is formalized")
    print("by large deviation theory (e.g., Sanov's theorem), which states that the probability of")
    print("observing a system in a non-equilibrium configuration is exponentially small.")
    print("For N bosons distributed in energy levels, with n_i particles in level i having g_i states,")
    print("the number of microstates is W = Π [ (n_i + g_i - 1)! / (n_i! * (g_i - 1)!) ].")
    print("Maximizing ln(W) is equivalent to minimizing the large deviation rate function.")
    print("-" * 70)

    # --- Step 2: Maximizing Entropy with Constraints (Cramer-Chernoff Principle) ---
    print("2. Maximizing Entropy (Connection to Cramer-Chernoff Theorem)")
    print("We maximize S subject to a fixed mean energy U = Σ n_i * ε_i.")
    print("For photons, particle number is not conserved, so the chemical potential is zero.")
    print("This constrained optimization, handled by Lagrange multipliers, is a physical application")
    print("of the principles behind the Cramer-Chernoff theorem, which deals with deviations of a sample mean.")
    print("The result of maximizing S under the energy constraint yields the Bose-Einstein distribution:")
    print("<n_i> / g_i = 1 / (exp(ε_i / (k_B * T)) - 1)")
    print("where <n_i> is the average number of photons in state i, T is the temperature, and k_B is Boltzmann's constant.")
    print("-" * 70)
    
    # --- Step 3: Calculating Mean Energy (U) ---
    print("3. Calculating the Mean Energy (U)")
    print("To find the total energy U, we integrate the energy per state over the density of states g(ε).")
    print("For photons (spin-1 bosons, 2 polarizations) in a volume V, g(ε) = (8 * π * V / (h*c)^3) * ε^2.")
    print("U = ∫ [0, ∞] ε * (<n(ε)>/g(ε)) * g(ε) dε")
    print("U = ∫ [0, ∞] ε * (1 / (exp(ε/(k_B*T)) - 1)) * ((8*π*V)/(h*c)^3 * ε^2) dε")
    print("This integral resolves to:")
    print("∫ [0, ∞] x^3 / (e^x - 1) dx = π^4 / 15, after substitution x = ε/(k_B*T).")
    print("This leads to the Stefan-Boltzmann law for total energy.")
    print("\n--- EQUILIBRIUM MEAN ENERGY (U) ---")
    print("The final expression for the mean energy U is:")
    print("U = (8 * π^5 * V * (k_B*T)^4) / (15 * (h*c)^3)")
    print("\nExplicitly, the numbers and symbols in the formula are:")
    print("U = (constant_1 * (pi**power_1) * V * (k_B**power_2) * (T**power_2)) / (constant_2 * (h**power_3) * (c**power_3))")
    print(f"  - constant_1 = {8}")
    print(f"  - power_1 (for pi) = {5}")
    print(f"  - constant_2 = {15}")
    print(f"  - power_2 (for k_B and T) = {4}")
    print(f"  - power_3 (for h and c) = {3}")
    print("  - V, k_B, T, h, c are Volume, Boltzmann constant, Temperature, Planck's constant, and speed of light.")
    print("-" * 70)
    
    # --- Step 4: Calculating Entropy (S) ---
    print("4. Calculating the Entropy (S)")
    print("The entropy can be found from the thermodynamic relation S = (U - F) / T,")
    print("where F is the Helmholtz free energy. F is calculated from the partition function Z:")
    print("F = -k_B * T * ln(Z) = -k_B * T * ∫ g(ε) * ln(1 - exp(-ε/(k_B*T))) dε")
    print("Solving this integral gives F = -U / 3.")
    print("Therefore, S = (U - (-U/3)) / T = (4/3) * U / T.")
    print("Substituting the expression for U:")
    
    print("\n--- EQUILIBRIUM ENTROPY (S) ---")
    print("The final expression for the entropy S is:")
    print("S = (32 * π^5 * V * k_B^4 * T^3) / (45 * (h*c)^3)")
    print("\nExplicitly, the numbers and symbols in the formula are:")
    print("S = (constant_1 * (pi**power_1) * V * (k_B**power_2) * (T**power_3)) / (constant_2 * (h**power_4) * (c**power_4))")
    print(f"  - constant_1 = {32}")
    print(f"  - power_1 (for pi) = {5}")
    print(f"  - constant_2 = {45}")
    print(f"  - power_2 (for k_B) = {4}")
    print(f"  - power_3 (for T) = {3}")
    print(f"  - power_4 (for h and c) = {3}")
    print("  - V, k_B, T, h, c are Volume, Boltzmann constant, Temperature, Planck's constant, and speed of light.")
    print("-" * 70)

if __name__ == '__main__':
    derive_bose_equilibrium()