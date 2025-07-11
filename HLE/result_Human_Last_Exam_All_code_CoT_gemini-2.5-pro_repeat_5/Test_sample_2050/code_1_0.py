import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given Hamiltonian.
    
    The problem describes N=7 coupled cavities in a ring with M=4 photons, in the limit of U -> infinity.
    This corresponds to a system of M=4 non-interacting fermions on a lattice of N=7 sites.
    The Hamiltonian for the fermions is H = M*omega + H_hop, where H_hop is the tight-binding Hamiltonian.
    The single-particle energies are given by epsilon_k = -2*J*cos(k), with k = 2*pi*m/N.
    For N=7, the integer m runs from -3 to 3.
    
    The energy levels in increasing order are for m = 0, (+-1), (+-2), (+-3).
    To find the ground state for M=4 fermions, we fill the 4 lowest energy states:
    - 1 fermion in the m=0 state.
    - 2 fermions in the m=+1 and m=-1 states.
    - 1 fermion in the m=+2 state (or m=-2, they are degenerate).
    
    The total ground state energy is the sum of these single-particle energies plus the constant on-site energy.
    E_g = M*omega + (epsilon_0 + epsilon_1 + epsilon_-1 + epsilon_2)
    """
    N = 7
    M = 4

    # Momenta are k_m = 2*pi*m/N. We need m=0, 1, -1, 2.
    # Note that cos(2*pi*(-1)/N) = cos(2*pi*1/N)
    k0 = 0
    k1 = 2 * np.pi * 1 / N
    k2 = 2 * np.pi * 2 / N
    
    # Calculate cosine terms
    cos_k0 = np.cos(k0)
    cos_k1 = np.cos(k1)
    cos_k2 = np.cos(k2)

    # The hopping energy is -2*J * (sum of cosines of occupied k-states)
    # E_hop = -2*J*(cos(k0) + cos(k1) + cos(k-1) + cos(k2))
    # E_hop = -2*J*(cos(k0) + 2*cos(k1) + cos(k2))
    coeff_J = -2 * (cos_k0 + 2 * cos_k1 + cos_k2)

    print("The ground state energy E_g is given by the formula:")
    print(f"E_g = {M}*omega - 2*J*(cos(2*pi*0/{N}) + 2*cos(2*pi*1/{N}) + cos(2*pi*2/{N}))")
    print("\nCalculating the numerical values:")
    print(f"cos(0) = {cos_k0:.5f}")
    print(f"cos(2*pi/{N}) = {cos_k1:.5f}")
    print(f"cos(4*pi/{N}) = {cos_k2:.5f}")
    
    print("\nSubstituting these values, the final expression for the ground state energy is:")
    print(f"E_g = {M} * omega + ({coeff_J:.5f}) * J")

calculate_ground_state_energy()