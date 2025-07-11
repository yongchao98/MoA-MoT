import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of hard-core bosons on a ring.

    The problem is mapped to non-interacting spinless fermions on a lattice. The
    ground state energy is the sum of the M lowest single-particle energies.
    """
    # System parameters
    N = 7  # Number of cavities (sites)
    M = 4  # Number of photons (fermions)

    # Calculate the single-particle energies (in units of J)
    # The energy levels are given by -2 * cos(2*pi*m/N)
    single_particle_energies_J = []
    for m in range(N):
        energy_coeff = -2 * np.cos(2 * np.pi * m / N)
        single_particle_energies_J.append(energy_coeff)

    # Sort the energies to find the lowest ones
    single_particle_energies_J.sort()

    # The ground state energy from the hopping term is the sum of the M lowest energies
    ground_state_hopping_energies_J = single_particle_energies_J[:M]
    total_hopping_energy_coeff_J = sum(ground_state_hopping_energies_J)

    # The total ground state energy is E = M*omega + E_hopping
    # We will print the final equation with all the numbers.

    print("The ground state energy E_gs is determined by filling the 4 lowest single-particle energy levels.")
    print("The total energy is E_gs = 4 * omega + E_hopping.")
    print("\nThe hopping energy is the sum of the four lowest single-particle energy contributions:")
    
    # Building the string for the equation
    equation_parts = [f"{e:.5f}" for e in ground_state_hopping_energies_J]
    equation_str = " + ".join(equation_parts)
    print(f"E_hopping = ({' + '.join(f'({e:.5f})' for e in ground_state_hopping_energies_J)}) * J")
    print(f"E_hopping = ({total_hopping_energy_coeff_J:.5f}) * J")

    print("\nThus, the final ground state energy is:")
    print(f"E_gs = {M} * omega + ({total_hopping_energy_coeff_J:.5f}) * J")

calculate_ground_state_energy()