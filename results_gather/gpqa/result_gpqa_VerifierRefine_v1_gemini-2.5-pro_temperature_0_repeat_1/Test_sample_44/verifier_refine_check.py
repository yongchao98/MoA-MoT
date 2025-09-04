import itertools

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the energy levels of four spin-1/2 particles
    in a 1D infinite potential well.
    """
    
    # --- Problem Definition ---
    num_particles = 4
    # Spin-1/2 particles are fermions, so max 2 particles per energy level (Pauli Exclusion Principle)
    max_occupancy_per_level = 2
    # The answer provided by the LLM to be checked.
    llm_answer = (10, 15, 18)

    # --- Physics Calculation ---
    
    # Single-particle energy for level n is n^2 * E. We calculate in units of E.
    def single_particle_energy(n):
        return n**2

    # We need to find all ways to place 4 particles into energy levels.
    # This is equivalent to finding integer partitions of 4, but with a twist:
    # we are partitioning particles into groups, and each group goes to a level.
    # A more direct approach is to generate all possible states.
    
    # We can generate states by considering all combinations of 4 levels chosen from a larger set,
    # allowing for replacement (since multiple particles can be in the same level).
    # We'll consider levels up to a reasonable limit (e.g., n=7) which is more than enough
    # to find the first few excited states.
    
    levels_to_consider = range(1, 8)
    all_possible_states = set()

    # Generate combinations with replacement to represent the levels of the 4 particles
    for state_levels in itertools.combinations_with_replacement(levels_to_consider, num_particles):
        # Check if the state is valid according to the Pauli Exclusion Principle
        is_valid = True
        level_counts = {}
        for level in state_levels:
            level_counts[level] = level_counts.get(level, 0) + 1
        
        for count in level_counts.values():
            if count > max_occupancy_per_level:
                is_valid = False
                break
        
        if is_valid:
            # Add the valid state (as a sorted tuple to ensure uniqueness) to our set
            all_possible_states.add(tuple(sorted(state_levels)))

    # Calculate the total energy for each unique valid state
    state_energies = set()
    for state in all_possible_states:
        total_energy = sum(single_particle_energy(n) for n in state)
        state_energies.add(total_energy)
        
    # Sort the unique energies to find the ground and excited states
    sorted_energies = sorted(list(state_energies))

    # The first three energies are the ground, first excited, and second excited states
    if len(sorted_energies) < 3:
        return "Error: Could not calculate at least three distinct energy states."
        
    calculated_energies = tuple(sorted_energies[:3])

    # --- Verification ---
    
    # Ground State Configuration: 2 particles in n=1, 2 in n=2. Energy = 2*(1^2) + 2*(2^2) = 10.
    # First Excited State Configuration: 2 in n=1, 1 in n=2, 1 in n=3. Energy = 2*(1^2) + 1*(2^2) + 1*(3^2) = 15.
    # Second Excited State Configuration: 1 in n=1, 2 in n=2, 1 in n=3. Energy = 1*(1^2) + 2*(2^2) + 1*(3^2) = 18.
    
    expected_energies = (10, 15, 18)

    if calculated_energies != expected_energies:
        return (f"Calculation mismatch. The code calculated the first three energy states to be {calculated_energies}E, "
                f"but the theoretical values are {expected_energies}E.")

    if calculated_energies == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The LLM's answer is {llm_answer}E, but the calculated correct answer is {calculated_energies}E.\n"
                f"Reasoning:\n"
                f"- Ground State: Particles fill the lowest available levels. Configuration: (2 particles in n=1, 2 in n=2). Energy = 2*1^2 + 2*2^2 = 10E.\n"
                f"- 1st Excited State: Promote one particle from the highest occupied level (n=2) to the lowest unoccupied level (n=3). Configuration: (2 in n=1, 1 in n=2, 1 in n=3). Energy = 2*1^2 + 1*2^2 + 1*3^2 = 15E.\n"
                f"- 2nd Excited State: The next lowest energy configuration. Comparing possibilities, promoting one particle from n=1 to n=3 gives the lowest energy increase after the first excitation. Configuration: (1 in n=1, 2 in n=2, 1 in n=3). Energy = 1*1^2 + 2*2^2 + 1*3^2 = 18E.")

# Run the check
result = check_correctness()
print(result)