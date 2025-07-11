def explain_bec_chemical_potential():
    """
    This function explains the fundamental limit on the chemical potential for bosons
    in the context of Bose-Einstein condensation (BEC).
    """

    print("Step-by-step derivation of the limit on chemical potential (μ) for bosons:")
    print("-" * 70)

    # Step 1: The Bose-Einstein distribution
    print("1. The starting point is the Bose-Einstein distribution for the number of particles n(ε) in a state with energy ε:")
    print("   n(ε) = 1 / (exp((ε - μ) / k_B*T) - 1)\n")

    # Step 2: The physical requirement
    print("2. For n(ε) to be a positive number, the denominator must be positive.")
    print("   This means: exp((ε - μ) / k_B*T) - 1 > 0")
    print("   Which simplifies to: ε - μ > 0, or μ < ε\n")

    # Step 3: The most restrictive limit
    print("3. This condition must hold for ALL energy states. The most restrictive case is for the lowest energy, the ground state (ε₀).")
    print("   Therefore, the fundamental limit on the chemical potential is:")
    print("   μ < ε₀\n")

    # Step 4: Condition for condensation
    print("4. BEC occurs when a macroscopic number of particles occupies the ground state (n(ε₀) -> ∞).")
    print("   This requires the denominator for n(ε₀) to approach zero: exp((ε₀ - μ) / k_B*T) - 1 -> 0")
    print("   This happens when μ gets infinitesimally close to ε₀.")
    print("   In the condensed phase (T ≤ T_c), we consider the chemical potential to be fixed at the ground state energy:")
    print("   Final Equation: μ = ε₀\n")

    # Step 5: Evaluating the correct option
    print("5. Let's analyze option C: 'μ must be equal to the chemical potential of a non-interacting Bose gas at zero temperature.'")
    print("   At T=0, all non-interacting bosons are in the ground state. The chemical potential is thus μ(T=0) = ε₀.")
    print("   This matches our conclusion from step 4. So, the condition for being in the BEC phase (μ = ε₀) is correctly described by option C.")
    print("-" * 70)

# Run the explanation
explain_bec_chemical_potential()