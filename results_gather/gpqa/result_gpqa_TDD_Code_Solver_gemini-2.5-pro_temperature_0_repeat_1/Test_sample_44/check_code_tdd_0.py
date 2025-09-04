import itertools

def check_answer():
    """
    This function checks the correctness of the given answer for the energies of
    four spin-1/2 particles in a 1D infinite potential well.

    The energy of a single particle in state n is E_n = n^2 * E.
    Spin-1/2 particles are fermions, so by the Pauli Exclusion Principle, a maximum
    of two particles (one spin-up, one spin-down) can occupy any given energy level n.

    The function calculates the ground, first excited, and second excited state energies
    and compares them to the provided answer B) 10E, 15E, 18E.
    """

    # The answer to check, from option B
    llm_answer = (10, 15, 18)

    # Helper function to calculate the total energy of a system configuration.
    # A configuration is a list of occupied quantum numbers n for the 4 particles.
    # e.g., [1, 1, 2, 2] for the ground state.
    def calculate_total_energy(config):
        # The energy is the sum of n^2 for each particle, in units of E.
        return sum(n**2 for n in config)

    # --- Ground State ---
    # To find the lowest energy, we fill the lowest available energy levels.
    # Level n=1 can hold 2 particles.
    # Level n=2 can hold the other 2 particles.
    ground_state_config = [1, 1, 2, 2]
    ground_energy = calculate_total_energy(ground_state_config)

    if ground_energy != llm_answer[0]:
        return f"Incorrect ground state energy. The lowest energy configuration is two particles in n=1 and two in n=2, giving a total energy of 2*(1^2) + 2*(2^2) = 10E. The answer provided {llm_answer[0]}E."

    # --- Excited States ---
    # To find the excited states, we need to find all possible valid configurations
    # and sort them by energy. A valid configuration has 4 particles and respects
    # the Pauli principle (max 2 particles per level).

    # We can generate configurations by considering partitions of 4 particles into levels.
    # We'll check levels up to n=5, which is sufficient for the first few excited states.
    possible_energies = set()
    max_n = 6 # A sufficiently large quantum number to find the first few states
    
    # Iterate through all combinations of 4 levels chosen from a list with duplicates
    # representing the available states (e.g., [1,1,2,2,3,3,...])
    available_states = []
    for i in range(1, max_n + 1):
        available_states.extend([i, i]) # Each level n is available for two fermions

    for config in itertools.combinations(available_states, 4):
        possible_energies.add(calculate_total_energy(config))

    # Sort the unique energies to find the ground, 1st, 2nd, etc. states
    sorted_energies = sorted(list(possible_energies))

    # The first element is the ground state energy
    calculated_ground = sorted_energies[0]
    # The second element is the first excited state energy
    calculated_first_excited = sorted_energies[1]
    # The third element is the second excited state energy
    calculated_second_excited = sorted_energies[2]

    calculated_answer = (calculated_ground, calculated_first_excited, calculated_second_excited)

    # Compare the calculated tuple with the LLM's answer tuple
    if calculated_answer != llm_answer:
        return (f"The calculated energies {calculated_answer}E do not match the provided answer {llm_answer}E. "
                f"Ground state should be {calculated_ground}E. "
                f"First excited state should be {calculated_first_excited}E. "
                f"Second excited state should be {calculated_second_excited}E.")

    return "Correct"

# Execute the check
result = check_answer()
print(result)