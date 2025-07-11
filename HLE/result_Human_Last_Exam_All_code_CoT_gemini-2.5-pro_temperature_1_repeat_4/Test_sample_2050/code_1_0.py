import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given Hamiltonian.
    """
    # System parameters
    N = 7
    num_photons = 4

    # The ground state energy E_gs has two components:
    # 1. A constant offset from the on-site energy term: num_photons * omega
    # 2. A contribution from the hopping term: E_hopping
    # E_gs = 4*omega + E_hopping

    # For U -> infinity, we have hard-core bosons, equivalent to non-interacting fermions.
    # The ground state is found by filling the lowest single-particle energy levels.
    # The single-particle energies are E_m = -2*J*cos(2*pi*m/N).

    # For N=7, the 4 lowest energy modes correspond to m=0, m=1, m=6, and m=2.
    # E_hopping = E_0 + E_1 + E_6 + E_2
    # E_hopping = -2*J * (cos(0) + cos(2*pi/7) + cos(12*pi/7) + cos(4*pi/7))
    # Since cos(12*pi/7) is the same as cos(2*pi/7), this simplifies to:
    # E_hopping = -2*J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))

    # Calculate the values of the cosine terms
    cos_0 = np.cos(0)
    cos_2pi_7 = np.cos(2 * np.pi / 7)
    cos_4pi_7 = np.cos(4 * np.pi / 7)

    # Calculate the term in the parenthesis
    sum_of_cosines = cos_0 + 2 * cos_2pi_7 + cos_4pi_7

    # Calculate the final coefficient for the J term
    J_coeff = -2 * sum_of_cosines

    # Print the step-by-step derivation of the final equation
    print("The ground state energy E_gs is given by the sum of the on-site energy and the hopping energy:")
    print(f"E_gs = {num_photons}*omega + E_hopping")
    
    print("\nThe hopping energy E_hopping is calculated by summing the 4 lowest single-particle energies:")
    print("E_hopping = -2*J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))")
    
    print("\nPlugging in the numerical values for the cosines:")
    print(f"cos(0) = {cos_0:.4f}")
    print(f"cos(2*pi/7) = {cos_2pi_7:.4f}")
    print(f"cos(4*pi/7) = {cos_4pi_7:.4f}")
    
    print("\nSubstituting these values into the expression for E_hopping:")
    print(f"E_hopping = -2*J * ({cos_0:.4f} + 2*{cos_2pi_7:.4f} + {cos_4pi_7:.4f})")
    print(f"E_hopping = -2*J * ({sum_of_cosines:.4f})")
    print(f"E_hopping = {J_coeff:.4f}*J")
    
    print("\nTherefore, the final equation for the ground state energy is:")
    print(f"E_gs = {num_photons} * omega + ({J_coeff:.4f}) * J")

calculate_ground_state_energy()