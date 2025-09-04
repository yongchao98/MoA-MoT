def check_correctness():
    """
    Checks the correctness of the calculated energies for four spin-1/2 particles
    in a 1D infinite potential well.

    The code calculates the energies for the ground state, first excited state,
    and second excited state based on first principles and compares them to the
    energies given in the proposed answer (10E, 15E, 18E).
    """

    # The energies from the LLM's answer, in units of E.
    llm_answer_energies = {
        "ground": 10,
        "first_excited": 15,
        "second_excited": 18
    }

    # Helper function to calculate the total energy of a system configuration.
    # A configuration is given as a list of occupation numbers for each level n=1, 2, 3...
    # e.g., [2, 2] means 2 particles in n=1, 2 particles in n=2.
    def calculate_total_energy(config):
        total_e = 0
        for i, num_particles in enumerate(config):
            n = i + 1  # Energy levels are n=1, 2, 3...
            # Single particle energy is n^2 * E
            total_e += num_particles * (n**2)
        return total_e

    # 1. Ground State Calculation
    # To minimize energy, fill the lowest levels first.
    # 4 particles, max 2 per level.
    # Level n=1: 2 particles
    # Level n=2: 2 particles
    ground_config = [2, 2]
    calculated_ground_energy = calculate_total_energy(ground_config)

    # Check ground state
    if calculated_ground_energy != llm_answer_energies["ground"]:
        return (f"Incorrect ground state energy. "
                f"The Pauli exclusion principle dictates filling the lowest energy levels first (2 particles in n=1, 2 in n=2). "
                f"This gives a calculated energy of 2*(1^2)E + 2*(2^2)E = {calculated_ground_energy}E, "
                f"but the answer provided is {llm_answer_energies['ground']}E.")

    # 2. Excited States Calculation
    # Excited states are found by promoting one or more particles to higher, unoccupied levels.
    # We list possible low-energy excitations from the ground state [2, 2].
    
    # Excitation A: Promote one particle from n=2 to n=3.
    # Config: [2, 1, 1] -> 2 in n=1, 1 in n=2, 1 in n=3
    energy_A = calculate_total_energy([2, 1, 1]) # 2*1 + 1*4 + 1*9 = 15

    # Excitation B: Promote one particle from n=1 to n=3.
    # Config: [1, 2, 1] -> 1 in n=1, 2 in n=2, 1 in n=3
    energy_B = calculate_total_energy([1, 2, 1]) # 1*1 + 2*4 + 1*9 = 18

    # Excitation C (two-particle): Promote both particles from n=2 to n=3.
    # Config: [2, 0, 2] -> 2 in n=1, 2 in n=3
    energy_C = calculate_total_energy([2, 0, 2]) # 2*1 + 2*9 = 20

    # Excitation D: Promote one particle from n=2 to n=4.
    # Config: [2, 1, 0, 1] -> 2 in n=1, 1 in n=2, 1 in n=4
    energy_D = calculate_total_energy([2, 1, 0, 1]) # 2*1 + 1*4 + 1*16 = 22

    # Collect all unique excited state energies and sort them to find the order.
    excited_energies = sorted(list(set([energy_A, energy_B, energy_C, energy_D])))

    # Check first excited state
    calculated_first_excited = excited_energies[0]
    if calculated_first_excited != llm_answer_energies["first_excited"]:
        return (f"Incorrect first excited state energy. "
                f"The lowest energy excitation from the ground state is {calculated_first_excited}E, "
                f"but the answer provided is {llm_answer_energies['first_excited']}E.")

    # Check second excited state
    calculated_second_excited = excited_energies[1]
    if calculated_second_excited != llm_answer_energies["second_excited"]:
        return (f"Incorrect second excited state energy. "
                f"The second-lowest energy excitation from the ground state is {calculated_second_excited}E, "
                f"but the answer provided is {llm_answer_energies['second_excited']}E.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)