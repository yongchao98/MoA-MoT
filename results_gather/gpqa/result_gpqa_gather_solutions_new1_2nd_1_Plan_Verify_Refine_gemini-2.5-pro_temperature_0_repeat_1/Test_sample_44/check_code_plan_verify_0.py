import collections

def check_correctness():
    """
    This function checks the correctness of the given answer for the quantum mechanics problem.

    The problem asks for the energies of the ground state, first excited state, and
    second excited state for four identical spin-1/2 particles in a 1D infinite potential well.

    Key Principles:
    1.  Single-particle energy levels: E_n = n^2 * E, where E is a constant.
    2.  Pauli Exclusion Principle: Since the particles are spin-1/2 (fermions), each
        energy level 'n' can hold a maximum of two particles.
    """

    # The final answer provided by the LLM is 'D'.
    llm_answer_key = 'D'

    # The options given in the question.
    options = {
        'A': (30, 39, 50),
        'B': (4, 10, 50),
        'C': (4, 10, 15),
        'D': (10, 15, 18)
    }

    # --- Step 1: Calculate the Ground State Energy ---
    # To find the ground state, we fill the lowest energy levels first, with two
    # particles per level.
    # Configuration: 2 particles in n=1, 2 particles in n=2.
    # Energy = (2 * 1^2 * E) + (2 * 2^2 * E) = (2 * 1 + 2 * 4) * E = (2 + 8) * E = 10E.
    e_ground = 2 * (1**2) + 2 * (2**2)

    # --- Step 2: Calculate the First Excited State Energy ---
    # This is the lowest energy excitation from the ground state. We promote one
    # particle from the highest occupied level (n=2) to the lowest unoccupied level (n=3).
    # Configuration: 2 particles in n=1, 1 particle in n=2, 1 particle in n=3.
    # Energy = (2 * 1^2 * E) + (1 * 2^2 * E) + (1 * 3^2 * E) = (2 + 4 + 9) * E = 15E.
    e_first_excited = 2 * (1**2) + 1 * (2**2) + 1 * (3**2)

    # --- Step 3: Calculate the Second Excited State Energy ---
    # This is the next lowest energy state. We must compare the energies of different
    # possible low-energy excitations from the ground state.
    
    # Possibility A (already calculated): Promote n=2 -> n=3. Total Energy = 15E.
    
    # Possibility B: Promote n=1 -> n=3.
    # Configuration: 1 particle in n=1, 2 in n=2, 1 in n=3.
    # Energy = (1 * 1^2) + (2 * 2^2) + (1 * 3^2) = 1 + 8 + 9 = 18E.
    energy_b = 1 * (1**2) + 2 * (2**2) + 1 * (3**2)

    # Possibility C: Promote both particles from n=2 -> n=3.
    # Configuration: 2 particles in n=1, 2 in n=3.
    # Energy = (2 * 1^2) + (2 * 3^2) = 2 + 18 = 20E.
    energy_c = 2 * (1**2) + 2 * (3**2)

    # Possibility D: Promote one particle from n=2 -> n=4.
    # Configuration: 2 particles in n=1, 1 in n=2, 1 in n=4.
    # Energy = (2 * 1^2) + (1 * 2^2) + (1 * 4^2) = 2 + 4 + 16 = 22E.
    energy_d = 2 * (1**2) + 1 * (2**2) + 1 * (4**2)

    # The list of possible total energies is [10, 15, 18, 20, 22, ...].
    # The ground state is 10E.
    # The first excited state is the next lowest, 15E.
    # The second excited state is the next lowest after that, 18E.
    e_second_excited = energy_b

    # --- Step 4: Verify the Answer ---
    calculated_energies = (e_ground, e_first_excited, e_second_excited)
    
    # Get the values corresponding to the LLM's chosen option.
    llm_answer_values = options.get(llm_answer_key)

    if llm_answer_values is None:
        return f"The provided answer key '{llm_answer_key}' is not a valid option (A, B, C, or D)."

    if calculated_energies == llm_answer_values:
        return "Correct"
    else:
        reason = (
            f"The calculated energies are (Ground: {calculated_energies[0]}E, "
            f"1st Excited: {calculated_energies[1]}E, "
            f"2nd Excited: {calculated_energies[2]}E).\n"
            f"This corresponds to the tuple {calculated_energies}.\n"
            f"The LLM's answer is '{llm_answer_key}', which corresponds to the tuple {llm_answer_values}.\n"
            f"The calculated values do not match the values in the selected option."
        )
        return reason

# Execute the check and print the result.
result = check_correctness()
print(result)