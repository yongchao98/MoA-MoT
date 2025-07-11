import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for a system of 4 hard-core bosons
    on a 7-site ring.
    """
    # System parameters
    N = 7  # Number of cavities/sites
    num_photons = 4  # Number of photons/fermions

    # We assume J=1 to calculate the energy coefficient. The final energy will be in units of J.
    J = 1.0

    print("Step 1: In the limit U -> infinity, the system is equivalent to non-interacting fermions.")
    print(f"We need to find the ground state energy for {num_photons} fermions on a ring with {N} sites.")
    print("-" * 30)

    print("Step 2: Calculate the single-particle energy levels.")
    print(f"The energies are given by E_m = -2*J*cos(2*pi*m/N) for m = 0..{N-1}.")
    
    # Calculate all single-particle energy levels
    energies = []
    for m in range(N):
        k = 2 * np.pi * m / N
        energy = -2 * J * np.cos(k)
        energies.append(energy)

    # Sort the energies to find the lowest ones
    energies.sort()

    print("-" * 30)
    print(f"Step 3: Fill the lowest {num_photons} energy levels to find the ground state energy.")
    
    # Select the lowest `num_photons` energies
    ground_state_levels = energies[:num_photons]

    # Calculate the total ground state energy
    ground_state_energy = sum(ground_state_levels)

    print("\nThe final ground state energy (relative to 4*omega) is the sum of the energies of the occupied states.")
    print("The individual energy levels being filled are (in units of J):")
    # Join the numbers into a single string for the equation
    equation_numbers = " + ".join([f"({e:.4f})" for e in ground_state_levels])
    print(f"E_gs / J = {equation_numbers}")
    
    # Display the final result
    print("\nWhich sums to:")
    print(f"E_gs / J = {ground_state_energy:.4f}")
    
    print("\nNote: The total Hamiltonian is H = 4*omega + H'. The value calculated is the ground state energy of H', denoted as E_gs.")


if __name__ == '__main__':
    calculate_ground_state_energy()
