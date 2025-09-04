import re
from itertools import combinations

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the quantum mechanics problem.

    The problem asks for the energies of the ground state, first excited state, and
    second excited state for four identical spin-1/2 particles (fermions) in a
    1D infinite potential well.

    Key principles:
    1. Single-particle energy levels are E_n = n^2 * E, where n=1, 2, 3...
    2. Pauli Exclusion Principle: Each energy level 'n' can hold at most 2 particles.
    """

    # The final answer provided by the LLM.
    llm_answer_text = "<<<C>>>"

    # The options given in the question. Energies are in units of E.
    options = {
        "A": [4, 10, 15],
        "B": [4, 10, 50],
        "C": [10, 15, 18],
        "D": [30, 39, 50]
    }

    # --- Step 1: Parse the LLM's chosen option ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure to parse the LLM's answer. It should be in the format <<<X>>>."
    
    chosen_option_key = match.group(1)
    chosen_option_energies = options.get(chosen_option_key)

    # --- Step 2: Calculate the correct energies from first principles ---

    def calculate_total_energy(config):
        """Calculates the total energy for a given particle configuration."""
        # config is a dictionary {level_n: num_particles}
        total_energy = 0
        for n, num_particles in config.items():
            total_energy += num_particles * (n**2)
        return total_energy

    # The ground state has the 4 fermions in the lowest possible energy levels.
    # {level_1: 2 particles, level_2: 2 particles}
    ground_state_config = {1: 2, 2: 2}
    ground_state_energy = calculate_total_energy(ground_state_config)

    # To find excited states, we generate configurations by promoting particles
    # from the ground state and find their energies.
    
    # We will store unique energies of valid configurations in a set.
    # A configuration is valid if it has 4 particles and at most 2 per level.
    all_energies = {ground_state_energy}

    # Generate possible excited configurations by moving 1 or 2 particles.
    # This is sufficient to find the first few excited states.
    
    # Single-particle promotions from the ground state {1:2, 2:2}
    for occupied_level in [1, 2]:
        for empty_level in range(3, 7): # Promote to levels 3, 4, 5, 6
            config = ground_state_config.copy()
            config[occupied_level] -= 1
            config[empty_level] = config.get(empty_level, 0) + 1
            all_energies.add(calculate_total_energy(config))

    # Two-particle promotions from the ground state {1:2, 2:2}
    # e.g., move both from n=2 to n=3 -> {1:2, 3:2}
    config_2p_A = {1: 2, 3: 2}
    all_energies.add(calculate_total_energy(config_2p_A))
    
    # e.g., move one from n=1 to n=3 and one from n=2 to n=3 -> {1:1, 2:1, 3:2}
    config_2p_B = {1: 1, 2: 1, 3: 2}
    all_energies.add(calculate_total_energy(config_2p_B))

    # Sort the unique energies to get the ordered states
    sorted_energies = sorted(list(all_energies))
    
    # The first three are the ground, 1st excited, and 2nd excited states.
    calculated_energies = sorted_energies[:3]

    # --- Step 3: Compare the calculated energies with the chosen option ---
    if calculated_energies == chosen_option_energies:
        return "Correct"
    else:
        reason = ""
        if calculated_energies[0] != chosen_option_energies[0]:
            reason += f"The ground state energy is incorrect. The correct value is {calculated_energies[0]}E, but the answer implies {chosen_option_energies[0]}E. "
        if calculated_energies[1] != chosen_option_energies[1]:
            reason += f"The first excited state energy is incorrect. The correct value is {calculated_energies[1]}E, but the answer implies {chosen_option_energies[1]}E. "
        if calculated_energies[2] != chosen_option_energies[2]:
            reason += f"The second excited state energy is incorrect. The correct value is {calculated_energies[2]}E, but the answer implies {chosen_option_energies[2]}E. "
        
        return f"Incorrect. {reason.strip()}"

# Run the check
result = check_correctness()
print(result)