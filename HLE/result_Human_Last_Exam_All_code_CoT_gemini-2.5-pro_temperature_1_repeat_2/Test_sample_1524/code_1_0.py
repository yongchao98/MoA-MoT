import sys

def explain_chemical_potential_limit():
    """
    Explains the fundamental limit on the chemical potential for bosons
    in the context of Bose-Einstein condensation.
    """
    print("This script explains the fundamental limit on the chemical potential (μ) for a system of bosons.")
    print("-" * 80)

    # Step 1: Introduce the Bose-Einstein Distribution
    print("Step 1: The starting point is the Bose-Einstein distribution.")
    print("This function gives the average number of bosons, n(ε), in a single-particle quantum state with energy ε:")
    print("n(ε) = 1 / (exp[(ε - μ) / (k_B * T)] - 1)")
    print("Here, μ is the chemical potential, T is temperature, and k_B is the Boltzmann constant.\n")

    # Step 2: State the Physical Constraint
    print("Step 2: A core physical principle is that the number of particles in any state cannot be negative.")
    print("Therefore, n(ε) must be greater than or equal to zero for all possible energies ε.")
    print("For this to be true, the denominator of the fraction must be positive:\nexp[(ε - μ) / (k_B * T)] - 1 > 0\n")

    # Step 3: Derive the Inequality for the Chemical Potential
    print("Step 3: We can solve this inequality for the chemical potential μ.")
    print("exp[(ε - μ) / (k_B * T)] > 1")
    print("Taking the natural logarithm of both sides (since ln(x) is an increasing function):")
    print("(ε - μ) / (k_B * T) > 0")
    print("Since T and k_B are positive, this simplifies to:")
    print("ε - μ > 0  --->  μ < ε\n")

    # Step 4: Identify the Most Restrictive Limit
    print("Step 4: The condition μ < ε must hold for every single-particle state available to the bosons.")
    print("The system has a lowest possible energy state, called the ground state, with energy ε₀.")
    print("To satisfy the condition for all states, μ must be less than the energy of the ground state.")
    print("This establishes the fundamental upper limit on the chemical potential.\n")

    # Step 5: The Final Relationship
    print("Step 5: As the system is cooled towards the condensation temperature, μ approaches ε₀ from below.")
    print("In the condensed state (T ≤ T_c), μ is effectively pinned at this limit (μ = ε₀) to allow for macroscopic occupation of the ground state.")
    print("At T=0, all non-interacting bosons are in the ground state, and their chemical potential is exactly ε₀.")
    print("Therefore, the limit is the chemical potential of a non-interacting Bose gas at T=0.\n")

    # Final Output
    print("-" * 80)
    print("The final relationship showing the fundamental limit on the chemical potential is:")

    # As requested, printing each component of the final equation: μ <= ε₀
    # Using Unicode for better looking symbols
    chemical_potential_symbol = "\u03BC"
    less_than_equal_symbol = "\u2264"
    ground_state_energy_symbol = "\u03B5\u2080"
    print(f"Symbol for Chemical Potential: {chemical_potential_symbol}")
    print(f"Symbol for 'is less than or equal to': {less_than_equal_symbol}")
    print(f"Symbol for Ground State Energy: {ground_state_energy_symbol}")
    print("\nFull Equation:")
    print(f"{chemical_potential_symbol} {less_than_equal_symbol} {ground_state_energy_symbol}")

explain_chemical_potential_limit()
<<<C>>>