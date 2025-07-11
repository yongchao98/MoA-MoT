import numpy as np

def solve_hamiltonian():
    """
    Calculates the ground state energy for the given Hamiltonian.
    """
    # Parameters of the system
    N = 7  # Number of cavities
    num_photons = 4  # Number of photons

    # --- Step 1: Simplify the Hamiltonian ---
    # In the limit U -> infinity, each cavity can hold at most one photon (n_i = 0 or 1).
    # This is the hard-core boson limit, which in 1D maps to non-interacting fermions.
    # The Hamiltonian becomes H = omega * sum(n_i) + H_hop.

    # --- Step 2: Separate the energy terms ---
    # The on-site energy is a constant offset: E_onsite = num_photons * omega.
    # The hopping Hamiltonian H_hop needs to be solved.
    # For 4 photons (an even number), the equivalent fermions have periodic boundary conditions.

    # --- Step 3: Calculate single-particle energies ---
    # The single-particle energies for a particle on an N-site ring are E_m = -2J * cos(2*pi*m/N).
    # We calculate these energies in units of J.
    
    energies_J = []
    for m in range(N):
        k = 2 * np.pi * m / N
        energy = -2 * np.cos(k)
        energies_J.append(energy)

    # Sort the energies to find the lowest ones
    energies_J.sort()

    # --- Step 4: Calculate the ground state energy ---
    # The ground state energy of the hopping part is the sum of the lowest `num_photons` energies.
    ground_state_hop_J = sum(energies_J[:num_photons])

    # --- Step 5: Print the detailed calculation ---
    print("The total ground state energy E_g is the sum of the on-site energy and the hopping energy.")
    print(f"The on-site energy is a constant offset: {num_photons}*omega.")
    
    print("\nThe hopping energy is the sum of the lowest 4 single-particle energies.")
    print("The single-particle energies E_m/J are given by -2*cos(2*pi*m/7).")
    
    # The 4 lowest energies correspond to filling the states with m=0, m=1, m=-1 (or 6), and m=2.
    e0_J = energies_J[0]
    e1_J = energies_J[1] # and energies_J[2]
    e2_J = energies_J[3]
    
    print(f"The four lowest single-particle energies (in units of J) are:")
    print(f"  E_0/J   = {e0_J:.5f}")
    print(f"  E_1/J   = {e1_J:.5f}")
    print(f"  E_-1/J  = {energies_J[2]:.5f}")
    print(f"  E_2/J   = {e2_J:.5f}")
        
    print(f"\nThe total hopping energy is the sum of these: ({e0_J:.5f} + {e1_J:.5f} + {energies_J[2]:.5f} + {e2_J:.5f})*J = {ground_state_hop_J:.5f}*J")
    
    print("\nTherefore, the final ground state energy is:")
    print(f"E_g = {num_photons}*omega + {ground_state_hop_J:.5f}*J")

if __name__ == '__main__':
    solve_hamiltonian()