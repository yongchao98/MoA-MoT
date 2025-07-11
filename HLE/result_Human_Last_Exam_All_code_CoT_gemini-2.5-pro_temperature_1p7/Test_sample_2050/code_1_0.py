import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given Hamiltonian in the U -> infinity limit.
    """
    N = 7  # Number of cavities
    M = 4  # Number of photons

    # In the U -> infinity limit, we have hard-core bosons, equivalent to spinless fermions.
    # The ground state is found by filling the M lowest single-particle energy levels.
    # The single-particle eigenenergies of the hopping Hamiltonian on an N-site ring are:
    # E_m = -2J * cos(2*pi*m/N)

    # Let's list the single-particle energy levels for N=7, characterized by mode number m.
    # The allowed values for m are 0, +/-1, +/-2, +/-3.
    # The energy levels are doubly degenerate for m != 0.
    
    # We need to pick the M=4 states with the lowest energies.
    # The energy E_m is minimized when cos(2*pi*m/N) is maximized.
    # cos(x) is largest for x=0. Then for x closest to 0.
    # Ordering of cos values: cos(0) > cos(2*pi/7) > cos(4*pi/7) > cos(6*pi/7).
    # This corresponds to energy levels E_0 < E_{pm 1} < E_{pm 2} < E_{pm 3}.

    # We fill the lowest energy states:
    # 1. One photon in the m=0 state.
    # 2. Two photons in the degenerate m=+1 and m=-1 states.
    # 3. One photon in one of the degenerate m=+2 or m=-2 states.
    
    # The sum of the energies of these four occupied states gives the hopping energy.
    
    # Wavevectors (k = 2*pi*m/N) for the four lowest energy states
    k_vals = [2*np.pi*0/N, 2*np.pi*1/N, 2*np.pi*(-1)/N, 2*np.pi*2/N]
    
    # Calculate the coefficients for J from each term -2*cos(k)
    coeffs_J = [-2 * np.cos(k) for k in k_vals]
    total_coeff_J = sum(coeffs_J)

    print("The ground state energy of the system is given by the sum of the on-site energy and the hopping energy:")
    print(f"E_ground = {M}*omega + E_hopping")
    print("\nThe hopping energy is the sum of the energies of the M=4 lowest-energy single-particle states:")
    print("E_hopping = E(m=0) + E(m=1) + E(m=-1) + E(m=2)")
    print("The energy for a single particle in mode m is E(m) = -2*J*cos(2*pi*m/N).")
    print("\nE_hopping = J * [ -2*cos(2*pi*0/7) - 2*cos(2*pi*1/7) - 2*cos(2*pi*(-1)/7) - 2*cos(2*pi*2/7) ]")
    
    print("\nCalculating the numerical values for the coefficients of J:")
    print(f"E_hopping = J * [ ({coeffs_J[0]:.5f}) + ({coeffs_J[1]:.5f}) + ({coeffs_J[2]:.5f}) + ({coeffs_J[3]:.5f}) ]")
    
    print(f"\nSumming the coefficients:")
    print(f"E_hopping = {total_coeff_J:.5f} * J")
    
    print("\nThus, the final equation for the ground state energy is:")
    print(f"E_ground = {M}*omega + ({total_coeff_J:.5f})*J")


if __name__ == '__main__':
    calculate_ground_state_energy()