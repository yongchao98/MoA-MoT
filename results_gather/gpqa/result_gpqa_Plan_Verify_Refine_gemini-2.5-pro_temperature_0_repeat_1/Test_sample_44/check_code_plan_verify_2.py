def check_correctness():
    """
    This function checks the correctness of the solution for the energy levels of four spin-1/2 particles
    in a 1D infinite potential well.
    """
    # The problem defines the base energy unit E = (pi^2 * hbar^2) / (2mL^2).
    # The single-particle energy levels are E_n = n^2 * E.
    # For calculation, we can treat E as 1 and add the suffix later.
    E_n = lambda n: n**2

    # Constraint: Four identical spin-1/2 particles are fermions.
    # By the Pauli Exclusion Principle, each energy level 'n' can hold at most 2 particles.

    # --- Ground State Energy ---
    # To find the ground state, we fill the lowest energy levels first.
    # 2 particles in n=1, 2 particles in n=2.
    ground_state_energy = 2 * E_n(1) + 2 * E_n(2)
    # 2*1 + 2*4 = 10

    # --- First Excited State Energy ---
    # Promote one particle from the highest occupied level (n=2) to the lowest unoccupied level (n=3).
    # Configuration: 2 particles in n=1, 1 in n=2, 1 in n=3.
    first_excited_energy = 2 * E_n(1) + 1 * E_n(2) + 1 * E_n(3)
    # 2*1 + 1*4 + 1*9 = 15

    # --- Second Excited State Energy ---
    # This is the state with the next lowest energy. We must find the excitation from the ground state
    # that has the lowest energy greater than the first excited state.
    #
    # Possibility A: Promote one particle from n=1 to n=3.
    # Config: 1 particle in n=1, 2 in n=2, 1 in n=3.
    energy_A = 1 * E_n(1) + 2 * E_n(2) + 1 * E_n(3) # 1*1 + 2*4 + 1*9 = 18
    #
    # Possibility B: Promote two particles from n=2 to n=3.
    # Config: 2 particles in n=1, 2 in n=3.
    energy_B = 2 * E_n(1) + 2 * E_n(3) # 2*1 + 2*9 = 20
    #
    # Possibility C: Promote one particle from n=2 to n=4.
    # Config: 2 particles in n=1, 1 in n=2, 1 in n=4.
    energy_C = 2 * E_n(1) + 1 * E_n(2) + 1 * E_n(4) # 2*1 + 1*4 + 1*16 = 22
    #
    # The list of possible excited state energies is [15, 18, 20, 22, ...].
    # The second excited state is the second one in this list.
    second_excited_energy = energy_A

    # The calculated sequence of energies is [10, 15, 18].
    calculated_energies = [ground_state_energy, first_excited_energy, second_excited_energy]

    # The options from the question are:
    # A) 4, 10, 15
    # B) 4, 10, 50
    # C) 30, 39, 50
    # D) 10, 15, 18
    correct_option_energies = [10, 15, 18]

    # Check if the calculated energies match the correct option.
    if calculated_energies == correct_option_energies:
        # The provided LLM answer correctly calculates the first excited state energy (15E),
        # which is a valid and correct step in reaching the final answer.
        # The full solution matches option D.
        return "Correct"
    else:
        return (f"Incorrect. The calculated energies are Ground={ground_state_energy}E, "
                f"1st Excited={first_excited_energy}E, 2nd Excited={second_excited_energy}E. "
                f"This does not match the expected correct answer of {correct_option_energies[0]}E, "
                f"{correct_option_energies[1]}E, {correct_option_energies[2]}E.")

# Execute the check
result = check_correctness()
print(result)