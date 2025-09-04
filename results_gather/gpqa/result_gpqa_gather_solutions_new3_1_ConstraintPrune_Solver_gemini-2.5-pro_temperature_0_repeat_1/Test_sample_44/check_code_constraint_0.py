def check_answer():
    """
    This function checks the correctness of the answer for the given quantum mechanics problem.
    It calculates the energies for the ground, first, and second excited states of four
    spin-1/2 particles in a 1D infinite potential well and compares them to the provided answer.
    """

    # The energy of the nth level is E_n = n^2 * E. We work with the coefficients of E.
    def calculate_total_energy(configuration):
        """Calculates the total energy for a given particle configuration."""
        return sum(n**2 for n in configuration)

    # --- Ground State ---
    # For fermions (spin-1/2 particles), we fill the lowest energy levels first,
    # with a maximum of two particles per level (Pauli Exclusion Principle).
    # Ground state configuration: 2 particles in n=1, 2 particles in n=2.
    ground_config = (1, 1, 2, 2)
    calculated_ground_energy = calculate_total_energy(ground_config)

    # --- Excited States ---
    # We need to find the configurations with the next lowest energies.
    # We list possible low-energy configurations and their total energies.
    
    # Ground state: {2 in n=1, 2 in n=2}, Energy = 10
    # 1st excited state (promote one n=2 -> n=3): {2 in n=1, 1 in n=2, 1 in n=3}
    first_excited_config = (1, 1, 2, 3)
    
    # Possible configurations for the second excited state:
    # A) Promote one n=1 -> n=3: {1 in n=1, 2 in n=2, 1 in n=3}
    second_excited_config_A = (1, 2, 2, 3)
    # B) Promote both n=2 -> n=3: {2 in n=1, 2 in n=3}
    second_excited_config_B = (1, 1, 3, 3)
    # C) Promote one n=2 -> n=4: {2 in n=1, 1 in n=2, 1 in n=4}
    second_excited_config_C = (1, 1, 2, 4)

    # Calculate energies for all configurations
    energies = {
        calculate_total_energy(ground_config),
        calculate_total_energy(first_excited_config),
        calculate_total_energy(second_excited_config_A),
        calculate_total_energy(second_excited_config_B),
        calculate_total_energy(second_excited_config_C)
    }

    # Sort the unique energies to find the ground, 1st, and 2nd excited states
    sorted_energies = sorted(list(energies))
    
    calculated_first_excited_energy = sorted_energies[1]
    calculated_second_excited_energy = sorted_energies[2]

    # The final answer provided is D, which corresponds to the energies 10E, 15E, 18E.
    expected_energies = [10, 15, 18]
    calculated_energies = [calculated_ground_energy, calculated_first_excited_energy, calculated_second_excited_energy]

    if calculated_energies == expected_energies:
        return "Correct"
    else:
        reason = "The calculated energies do not match the expected answer D (10E, 15E, 18E).\n"
        if calculated_ground_energy != expected_energies[0]:
            reason += f"Calculated ground state energy is {calculated_ground_energy}E, but should be {expected_energies[0]}E. This is derived from the configuration {ground_config}.\n"
        if calculated_first_excited_energy != expected_energies[1]:
            reason += f"Calculated first excited state energy is {calculated_first_excited_energy}E, but should be {expected_energies[1]}E. This is derived from the configuration {first_excited_config}.\n"
        if calculated_second_excited_energy != expected_energies[2]:
            reason += f"Calculated second excited state energy is {calculated_second_excited_energy}E, but should be {expected_energies[2]}E. This is derived from the configuration {second_excited_config_A}, which has lower energy than other possibilities like {second_excited_config_B} (Energy={calculate_total_energy(second_excited_config_B)}E).\n"
        return f"Incorrect. {reason.strip()}"

# print(check_answer())