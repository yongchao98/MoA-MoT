def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons in a
    Bose-Einstein condensate.
    """
    # Define variables as strings for the explanation
    mu = "μ"
    epsilon = "ε"
    epsilon_0 = "ε_0"
    k_B_T = "k_B*T"
    n_epsilon = "<n_ε>"

    print("Step 1: The starting point is the Bose-Einstein distribution.")
    print(f"The average occupation number ({n_epsilon}) of a single-particle state with energy {epsilon} is given by:")
    print(f"  {n_epsilon} = 1 / (exp(({epsilon} - {mu}) / ({k_B_T})) - 1)\n")

    print("Step 2: Apply the physical constraint.")
    print("The occupation number must be a non-negative quantity for any state.")
    print(f"  {n_epsilon} >= 0")
    print("This implies that the denominator must be positive and real.\n")

    print("Step 3: Derive the mathematical condition on the chemical potential.")
    print("For the denominator to be positive:")
    print(f"  exp(({epsilon} - {mu}) / ({k_B_T})) - 1 > 0")
    print("  exp(({epsilon} - {mu}) / ({k_B_T})) > 1")
    print("Taking the natural log of both sides:")
    print(f"  ({epsilon} - {mu}) / ({k_B_T}) > 0")
    print("Since thermal energy is positive, this simplifies to:")
    print(f"  {epsilon} - {mu} > 0  =>  {mu} < {epsilon}\n")

    print("Step 4: Identify the most stringent limit.")
    print(f"The condition {mu} < {epsilon} must be true for all energy states available to the bosons.")
    print(f"The most restrictive condition is therefore set by the lowest possible energy state, the ground state ({epsilon_0}).")
    print(f"So, the fundamental limit is: {mu} <= {epsilon_0}\n")
    
    print("Step 5: Relate the limit to Bose-Einstein Condensation (BEC).")
    print("In the condensed phase (at temperatures T <= T_c), a macroscopic number of particles occupies the ground state.")
    print(f"For the occupation of the ground state ({epsilon_0}) to become very large, its denominator must approach zero.")
    print(f"This means {mu} must approach {epsilon_0} from below. In the presence of a condensate, the chemical potential is effectively 'pinned' at the ground state energy.\n")
    
    print("The final condition for the chemical potential in a system with a Bose-Einstein condensate is:")
    print(f"  {mu} = {epsilon_0}")
    print("\nThis value, the ground state energy, is also the chemical potential of a non-interacting Bose gas at T=0.")
    print("This directly corresponds to choice C.")

explain_chemical_potential_limit()
<<<C>>>