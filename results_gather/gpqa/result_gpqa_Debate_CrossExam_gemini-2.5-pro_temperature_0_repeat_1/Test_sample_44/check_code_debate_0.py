import collections

def check_answer():
    """
    Checks the correctness of the calculated energies for a system of four spin-1/2 particles
    in a 1D infinite potential well.

    The problem asks for the energies of the ground state, first excited state, and
    second excited state in units of E = (pi^2 * hbar^2) / (2mL^2).

    The single-particle energy levels are E_n = n^2 * E.
    Since the particles are spin-1/2 fermions, each energy level 'n' can hold a maximum
    of two particles (Pauli Exclusion Principle).
    """

    # The answer provided by the LLM, corresponding to option A.
    # Energies are in units of E.
    llm_answer = [10, 15, 18]

    # --- Calculation ---

    # Define a helper function for single-particle energy in units of E
    def single_particle_energy(n):
        return n**2

    # 1. Ground State (E_gs)
    # We fill the lowest available energy levels with the four particles.
    # Configuration: Two particles in n=1, two particles in n=2.
    e_gs = 2 * single_particle_energy(1) + 2 * single_particle_energy(2)
    # e_gs = 2 * 1 + 2 * 4 = 10

    # 2. First Excited State (E_1st)
    # To get the first excited state, we find the lowest energy configuration
    # above the ground state. This is achieved by promoting one particle from the
    # highest occupied level (n=2) to the lowest unoccupied level (n=3).
    # Configuration: Two in n=1, one in n=2, one in n=3.
    e_1st = 2 * single_particle_energy(1) + 1 * single_particle_energy(2) + 1 * single_particle_energy(3)
    # e_1st = 2 * 1 + 1 * 4 + 1 * 9 = 15

    # 3. Second Excited State (E_2nd)
    # We must find the next possible configuration with the lowest energy
    # above the first excited state. We need to consider the next possible promotions
    # from the ground state configuration (1,1,2,2).
    
    # Possibility A: Promote one particle from n=1 to n=3.
    # Config: (1, 2, 2, 3)
    energy_A = 1 * single_particle_energy(1) + 2 * single_particle_energy(2) + 1 * single_particle_energy(3)
    # energy_A = 1 * 1 + 2 * 4 + 1 * 9 = 18

    # Possibility B: Promote both particles from n=2 to n=3.
    # Config: (1, 1, 3, 3)
    energy_B = 2 * single_particle_energy(1) + 2 * single_particle_energy(3)
    # energy_B = 2 * 1 + 2 * 9 = 20

    # Possibility C: Promote one particle from n=2 to n=4.
    # Config: (1, 1, 2, 4)
    energy_C = 2 * single_particle_energy(1) + 1 * single_particle_energy(2) + 1 * single_particle_energy(4)
    # energy_C = 2 * 1 + 1 * 4 + 1 * 16 = 22

    # The second excited state has the lowest energy among all excited states
    # that is greater than the first excited state's energy.
    # The energies of the first few excited states are 15E, 18E, 20E, 22E...
    # So the second excited state energy is the minimum of these possibilities.
    e_2nd = min(energy_A, energy_B, energy_C)

    calculated_energies = [e_gs, e_1st, e_2nd]

    # --- Verification ---
    if calculated_energies == llm_answer:
        return "Correct"
    else:
        reason = f"The calculated energies {calculated_energies}E do not match the LLM's answer {llm_answer}E.\n"
        if e_gs != llm_answer[0]:
            reason += f"Ground state energy calculation is incorrect. Expected {llm_answer[0]}E, but calculated {e_gs}E.\n"
        if e_1st != llm_answer[1]:
            reason += f"First excited state energy calculation is incorrect. Expected {llm_answer[1]}E, but calculated {e_1st}E.\n"
        if e_2nd != llm_answer[2]:
            reason += f"Second excited state energy calculation is incorrect. Expected {llm_answer[2]}E, but calculated {e_2nd}E.\n"
            reason += f"The possible next excited state energies were {energy_A}E, {energy_B}E, and {energy_C}E. The minimum of these is {e_2nd}E."
        return reason

# Run the check
result = check_answer()
print(result)