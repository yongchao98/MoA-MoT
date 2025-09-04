import itertools
from collections import Counter

def check_fermion_well_energies():
    """
    Calculates the energies for the ground state, first excited state, and
    second excited state for four spin-1/2 particles in a 1D infinite well.
    """
    
    # The LLM's answer is option A, which corresponds to energies 10E, 15E, 18E.
    llm_answer_energies = [10, 15, 18]

    # --- Step 1: Systematically find the lowest energy configurations ---
    # We need to find combinations of 4 quantum numbers {n1, n2, n3, n4}
    # that give the lowest total energies.
    # Constraint: Any 'n' value can appear at most twice (Pauli Exclusion Principle).

    # We can generate all valid low-energy configurations and sort them.
    # Generating combinations of quantum numbers from a reasonable range (e.g., 1 to 6)
    # is sufficient to find the first few excited states.
    possible_energies = set()
    for config in itertools.combinations_with_replacement(range(1, 7), 4):
        # Check if the configuration is valid according to the Pauli principle.
        counts = Counter(config)
        if all(count <= 2 for count in counts.values()):
            # If valid, calculate the total energy in units of E.
            # Total Energy = (n1^2 + n2^2 + n3^2 + n4^2) * E
            energy = sum(n**2 for n in config)
            possible_energies.add(energy)

    # Sort the unique energies to find the lowest ones.
    sorted_energies = sorted(list(possible_energies))
    
    # The ground, first, and second excited state energies are the first three in the list.
    if len(sorted_energies) < 3:
        return "Error: Could not determine the first three energy states."
        
    calculated_energies = sorted_energies[:3]

    # --- Step 2: Manually verify the configurations for clarity ---
    # Ground State (E_gs): Fill the lowest available levels.
    # Two particles in n=1, two particles in n=2. Config: {1, 1, 2, 2}
    # E_gs = (1^2 + 1^2 + 2^2 + 2^2)E = (1 + 1 + 4 + 4)E = 10E
    gs_config = [1, 1, 2, 2]
    gs_energy = sum(n**2 for n in gs_config)

    # First Excited State (E_1st): Promote one particle with the smallest energy jump.
    # Promote one particle from n=2 to n=3. Config: {1, 1, 2, 3}
    # E_1st = (1^2 + 1^2 + 2^2 + 3^2)E = (1 + 1 + 4 + 9)E = 15E
    ex1_config = [1, 1, 2, 3]
    ex1_energy = sum(n**2 for n in ex1_config)

    # Second Excited State (E_2nd): Find the next lowest energy configuration.
    # We compare the next possible promotions from the ground state:
    #   - Promote n=1 to n=3 -> {1, 2, 2, 3}: E = (1+4+4+9)E = 18E
    #   - Promote n=2 to n=4 -> {1, 1, 2, 4}: E = (1+1+4+16)E = 22E
    #   - Promote both n=2 to n=3 -> {1, 1, 3, 3}: E = (1+1+9+9)E = 20E
    # The lowest of these is 18E. Config: {1, 2, 2, 3}
    ex2_config = [1, 2, 2, 3]
    ex2_energy = sum(n**2 for n in ex2_config)
    
    manual_energies = [gs_energy, ex1_energy, ex2_energy]

    # --- Step 3: Compare calculated results with the LLM's answer ---
    # First, ensure our systematic and manual calculations agree.
    if calculated_energies != manual_energies:
        return (f"Internal check failed. Systematic calculation {calculated_energies} "
                f"does not match manual verification {manual_energies}.")

    # Now, check if the calculated energies match the LLM's answer.
    if calculated_energies == llm_answer_energies:
        return "Correct"
    else:
        return (f"Incorrect. The calculated energies are {calculated_energies[0]}E, "
                f"{calculated_energies[1]}E, and {calculated_energies[2]}E. "
                f"The LLM's answer corresponds to {llm_answer_energies[0]}E, "
                f"{llm_answer_energies[1]}E, and {llm_answer_energies[2]}E.\n"
                f"Reasoning:\n"
                f"- Ground State: Configuration {{1,1,2,2}}, Energy = (1+1+4+4)E = 10E.\n"
                f"- 1st Excited State: Configuration {{1,1,2,3}}, Energy = (1+1+4+9)E = 15E.\n"
                f"- 2nd Excited State: Configuration {{1,2,2,3}}, Energy = (1+4+4+9)E = 18E.")

# Run the check
result = check_fermion_well_energies()
print(result)