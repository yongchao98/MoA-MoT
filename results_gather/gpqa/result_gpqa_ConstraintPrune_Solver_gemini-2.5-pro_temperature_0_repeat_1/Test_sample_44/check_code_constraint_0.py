import itertools

def check_particle_well_answer():
    """
    Checks the correctness of the answer for the energy states of four spin-1/2
    particles in a 1D infinite potential well.
    """
    # The provided answer D corresponds to energies 10E, 15E, 18E.
    # We will check if these are the correct ground, 1st, and 2nd excited state energies.
    answer_energies = [10, 15, 18]

    # --- Model the Physics ---
    NUM_PARTICLES = 4
    # For spin-1/2 fermions, the Pauli Exclusion Principle allows a maximum of
    # two particles (one spin-up, one spin-down) per energy level 'n'.
    MAX_OCCUPANCY_PER_LEVEL = 2

    # --- Systematically Find All Low-Lying Energy States ---

    # We generate all possible ways to place the 4 particles into the first few
    # energy levels. Checking up to n=6 is sufficient to find the lowest states.
    # A configuration is a tuple of the principal quantum numbers for the 4 particles.
    # e.g., (1, 1, 2, 2) means 2 particles in n=1 and 2 in n=2.
    max_n_to_check = 6
    possible_levels = range(1, max_n_to_check + 1)
    
    # Use itertools to generate all combinations with replacement.
    all_configs = itertools.combinations_with_replacement(possible_levels, NUM_PARTICLES)

    valid_energies = set()

    for config in all_configs:
        # Check if the configuration is valid according to the Pauli Exclusion Principle.
        is_valid = True
        # We only need to check the levels present in the current configuration.
        for level in set(config):
            if config.count(level) > MAX_OCCUPANCY_PER_LEVEL:
                is_valid = False
                break
        
        if is_valid:
            # If valid, calculate the total energy in units of E.
            # The energy for a particle in level n is n^2 * E.
            total_energy = sum(n**2 for n in config)
            valid_energies.add(total_energy)

    # Sort the unique energies to find the ground and excited states in order.
    calculated_energies_sorted = sorted(list(valid_energies))

    # We need at least three energy levels to check the answer.
    if len(calculated_energies_sorted) < 3:
        return (f"Error: The code could not find at least 3 unique energy states. "
                f"Found: {calculated_energies_sorted}. This indicates a problem with the checking logic.")

    # The first three energies are the ground, 1st excited, and 2nd excited states.
    ground_E = calculated_energies_sorted[0]
    first_excited_E = calculated_energies_sorted[1]
    second_excited_E = calculated_energies_sorted[2]
    
    calculated_answer = [ground_E, first_excited_E, second_excited_E]

    # --- Compare with the provided answer ---
    if calculated_answer == answer_energies:
        return "Correct"
    else:
        # Construct a detailed reason for the failure.
        reason = (
            f"The answer is incorrect. The calculated energies for the ground, first, and second excited states are "
            f"{calculated_answer[0]}E, {calculated_answer[1]}E, and {calculated_answer[2]}E respectively.\n"
            f"The provided answer corresponds to {answer_energies[0]}E, {answer_energies[1]}E, and {answer_energies[2]}E.\n\n"
            f"Detailed breakdown of correct states:\n"
            f"- Ground State: Lowest energy configuration is placing 2 particles in n=1 and 2 in n=2. "
            f"Configuration: (1,1,2,2). Energy = 2*(1^2) + 2*(2^2) = 2 + 8 = 10E.\n"
            f"- 1st Excited State: The next lowest energy is achieved by promoting one particle from n=2 to n=3. "
            f"Configuration: (1,1,2,3). Energy = 2*(1^2) + 1*(2^2) + 1*(3^2) = 2 + 4 + 9 = 15E.\n"
            f"- 2nd Excited State: The third lowest energy configuration is promoting one particle from n=1 to n=3 (from the ground state). "
            f"Configuration: (1,2,2,3). Energy = 1*(1^2) + 2*(2^2) + 1*(3^2) = 1 + 8 + 9 = 18E."
        )
        return reason

# Execute the check and print the result.
result = check_particle_well_answer()
print(result)