def explain_bec_chemical_potential():
    """
    Derives and explains the fundamental limit on the chemical potential (μ)
    for bosons undergoing Bose-Einstein condensation.
    """
    print("Step 1: The Bose-Einstein Distribution")
    print("The average number of bosons, n(ε), in a state with energy ε is given by the equation:")
    # The prompt requires printing numbers/symbols in the final equation.
    # We will print the symbolic equation.
    print("n(ε) = 1 / (exp((ε - μ) / (k_B * T)) - 1)")
    print("where μ is the chemical potential.\n")

    print("Step 2: The Physical Constraint")
    print("A key physical requirement is that the number of particles, n(ε), cannot be negative.")
    print("For n(ε) to be positive, its denominator must also be positive:")
    print("exp((ε - μ) / (k_B * T)) - 1 > 0\n")

    print("Step 3: Deriving the Condition for μ")
    print("The inequality simplifies to:")
    print("exp((ε - μ) / (k_B * T)) > 1")
    print("Taking the natural logarithm of both sides, we get:")
    print("(ε - μ) / (k_B * T) > 0")
    print("Since temperature T and the Boltzmann constant k_B are positive, this implies:")
    print("ε - μ > 0, which means μ < ε\n")

    print("Step 4: The Fundamental Limit")
    print("This condition, μ < ε, must hold true for all possible energy states ε.")
    print("Therefore, the chemical potential μ must be less than the lowest possible energy, the ground state energy, ε_0.")
    print("The ultimate limit that the chemical potential can approach is the ground state energy itself.")
    print("So, the fundamental limit is: μ ≤ ε_0\n")

    print("Step 5: Relation to Condensation")
    print("Bose-Einstein condensation is the macroscopic occupation of the ground state.")
    print("For n(ε_0) to become very large, the denominator for the ground state must approach zero.")
    print("This happens when μ gets infinitesimally close to ε_0 from below.\n")

    print("Step 6: Conclusion and Matching the Answer")
    print("The fundamental limit for μ is the ground state energy, ε_0.")
    print("By definition, the chemical potential of a non-interacting Bose gas at absolute zero temperature (T=0) is exactly equal to its ground state energy, as all particles occupy this state.")
    print("Thus, the limit μ ≤ ε_0 is correctly described by choice C.\n")


explain_bec_chemical_potential()
<<<C>>>