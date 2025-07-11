import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given problem.

    The problem of hard-core bosons on a ring maps to non-interacting fermions.
    The ground state energy is found by filling the lowest available single-particle
    energy levels.
    """
    # Number of cavities
    N = 7
    # Number of photons
    M = 4

    # The Hamiltonian for hard-core bosons is equivalent to that of free fermions.
    # For an even number of particles (M=4), periodic boundary conditions apply.
    # The single-particle energies are E_q = -2*J*cos(2*pi*q/N).

    # We need to fill the M=4 lowest energy states. The energy is lowest when q
    # is close to 0. The occupied states are q=0, q=+/-1, and one of q=+/-2.
    occupied_q = [0, 1, -1, 2]

    # The total energy is E_ground = M*omega + E_hopping.
    # The hopping energy is the sum of the single-particle energies of the
    # occupied states.
    # E_hop = sum over q_occ of E_q = -2*J * sum(cos(2*pi*q/N))
    
    sum_cos = 0
    for q in occupied_q:
        sum_cos += np.cos(2 * np.pi * q / N)
    
    # The total coefficient for the J term in the energy
    coeff_J = -2 * sum_cos

    # The problem asks to output the final equation with each number.
    # The final equation is E = M * omega + coeff_J * J
    
    print("The final equation for the ground state energy is of the form: E = A * omega + B * J")
    print("\nEach number in the final equation is:")
    print(f"A = {M}")
    print(f"B = {coeff_J:.5f}")

    print("\nThus, the ground state energy is:")
    print(f"E_ground = {M} * omega + ({coeff_J:.5f}) * J")

calculate_ground_state_energy()