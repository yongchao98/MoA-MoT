import sys

def analyze_chemical_potential():
    """
    Analyzes the fundamental limit on the chemical potential for bosons
    in Bose-Einstein condensation.
    """
    # Step 1: The Bose-Einstein Distribution
    print("Step 1: The starting point is the Bose-Einstein distribution function.")
    print("The average number of bosons, n(ε), in a quantum state with energy ε is:")
    print("n(ε) = 1 / (exp((ε - μ) / (k_B * T)) - 1)")
    print("where μ is the chemical potential, T is the temperature, and k_B is Boltzmann's constant.\n")

    # Step 2: The Physical Constraint
    print("Step 2: A core physical constraint is that the number of particles, n(ε), must be positive.")
    print("For n(ε) > 0, the denominator of the fraction must also be positive:")
    print("exp((ε - μ) / (k_B * T)) - 1 > 0\n")

    # Step 3: Deriving the Limit on μ
    print("Step 3: From this inequality, we can derive the limit on μ.")
    print("exp((ε - μ) / (k_B * T)) > 1")
    print("Taking the natural logarithm of both sides gives:")
    print("(ε - μ) / (k_B * T) > 0")
    print("Since temperature T and k_B are positive, this means:")
    print("ε - μ > 0  or  μ < ε")
    print("This must be true for ALL energy states ε. The most restrictive condition comes from the lowest possible energy state, the ground state energy, ε_0.")
    
    # Define the variables for the final equation symbolically
    chemical_potential = "μ"
    ground_state_energy = "ε_0"
    
    # The final equation for the fundamental limit
    print("\n--- The Fundamental Limit ---")
    print(f"The chemical potential must be less than the ground state energy.")
    # The following line prints the "equation" as requested
    print(f"Final Equation of the Limit: {chemical_potential} < {ground_state_energy}")
    print("---------------------------\n")

    # Step 4: Relation to Bose-Einstein Condensation
    print("Step 4: Understanding the chemical potential *in* the condensate.")
    print("BEC is the macroscopic occupation of the ground state (n(ε_0) -> ∞).")
    print("For n(ε_0) to be very large, its denominator must approach zero.")
    print("This means (ε_0 - μ) must approach zero from the positive side.")
    print(f"Therefore, as the system condenses (T <= T_c), the chemical potential gets 'pinned' and approaches the ground state energy from below: {chemical_potential} → {ground_state_energy}")
    print(f"In the idealized case of a non-interacting gas in the thermodynamic limit, for T ≤ T_c, we have exactly: {chemical_potential} = {ground_state_energy}\n")
    
    # Step 5: Evaluating the Answer Choices
    print("Step 5: Evaluating the given choices.")
    print("The question asks for the limit on μ *in* Bose-Einstein condensation, which means we are considering the situation at or below the critical temperature.")
    print("Let's analyze option C: 'μ must be equal to the chemical potential of a non-interacting Bose gas at zero temperature.'")
    print("The chemical potential of a non-interacting Bose gas at T=0 is, by definition, the ground state energy, ε_0.")
    print("Therefore, option C is a different way of stating the condition that μ = ε_0, which is the correct description for an ideal Bose gas in the condensed phase.")

analyze_chemical_potential()