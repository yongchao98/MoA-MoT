import math

def display_photon_gas_properties():
    """
    This function uses established formulas from statistical mechanics,
    which are derived using principles underpinned by large deviation theory,
    to describe the equilibrium state of a photon gas (Bose case of light quanta).
    """

    print("The equilibrium state of a photon gas is found by maximizing its entropy.")
    print("This leads to the following key results:\n")

    # --- 1. The Equilibrium Distribution (Planck Distribution) ---
    print("1. Equilibrium Distribution of Quanta:")
    print("The average number of photons <n(ε)> in a state with energy ε = ħω is given by the Planck distribution:")
    print("   <n(ε)> = 1 / (exp(ε / (k_B * T)) - 1)")
    print("Where:")
    print("   ε = energy of the quantum")
    print("   k_B = Boltzmann constant")
    print("   T = Temperature")
    print("-" * 50)

    # --- 2. Equilibrium Mean Energy <E> ---
    print("2. Equilibrium Mean Energy <E>:")
    print("Integrating ε * <n(ε)> over the density of states for photons in a volume V gives the total mean energy <E>.")
    print("The result is the Stefan-Boltzmann law for total energy:")
    print("\n   <E> = (π² * V * k_B⁴ / (15 * ħ³ * c³)) * T⁴")
    print("\nBreaking down the numerical and symbolic parts of the equation for <E>:")
    print(f"   - Numerical Factor: (π² / 15) ≈ {math.pi**2 / 15:.4f}")
    print("   - System/Fundamental Constants: V * (k_B⁴ / (ħ³ * c³))")
    print("   - Temperature Dependence: T⁴")
    print("\nWhere:")
    print("   π = Pi")
    print("   V = Volume of the container")
    print("   k_B = Boltzmann constant")
    print("   ħ = Reduced Planck constant")
    print("   c = Speed of light")
    print("   T = Temperature")
    print("-" * 50)

    # --- 3. Equilibrium Entropy S ---
    print("3. Equilibrium Entropy S:")
    print("Using thermodynamic relations, the entropy S of the photon gas at equilibrium is found to be:")
    print("\n   S = (4 * π² * V * k_B⁴ / (45 * ħ³ * c³)) * T³")
    print("\nBreaking down the numerical and symbolic parts of the equation for S:")
    print(f"   - Numerical Factor: (4π² / 45) ≈ {4 * math.pi**2 / 45:.4f}")
    print("   - System/Fundamental Constants: V * (k_B⁴ / (ħ³ * c³))")
    print("   - Temperature Dependence: T³")
    print("-" * 50)

    # --- 4. Relationship between E and S ---
    print("4. A simple and elegant relationship exists between the equilibrium energy and entropy:")
    print("\n   S = (4/3) * <E> / T")
    print("\nThis shows that entropy is directly proportional to the mean energy and inversely proportional to the temperature.")

# Execute the function to print the results
display_photon_gas_properties()

<<<S = (4/3) * <E> / T>>>