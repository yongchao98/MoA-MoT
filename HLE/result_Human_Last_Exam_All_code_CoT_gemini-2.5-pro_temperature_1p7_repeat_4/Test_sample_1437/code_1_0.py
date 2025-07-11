import numpy as np

def solve_lindhard_function():
    """
    Determines the numerical value of the Lindhard polarization function at q=0, omega=0.

    The Lindhard function Pi(q, omega) in the limit of zero momentum transfer (q -> 0)
    and zero frequency (omega = 0) is given by the negative of the density of states
    at the Fermi energy, g(epsilon_F). This includes the factor of 2 for spin.

    So, the key physical relationship is:
    Pi(q=0, omega=0) = -g(epsilon_F)

    The quantity g(epsilon_F) depends on material properties (specifically, the electron
    density), and it has units of [Energy]^-1 * [Volume]^-1. Therefore, Pi(0,0) itself
    is not a universal numerical constant.

    A universal, dimensionless numerical value is found by normalizing Pi(0,0) by the
    density of states g(epsilon_F). This ratio is constant for any 3D homogeneous
    electron gas.
    """

    # Define the parameters of the limit we are evaluating
    momentum_transfer_q = 0
    external_frequency_omega = 0

    # The equation for the normalized value is:
    # Normalized_Pi = Pi(q=0, omega=0) / g(epsilon_F)
    # Substituting Pi(q=0, omega=0) = -g(epsilon_F), we get:
    # Normalized_Pi = -g(epsilon_F) / g(epsilon_F)
    normalized_pi_value = -1.0

    print("--- Lindhard Polarization Function Evaluation ---")
    print("\nWe are calculating the Lindhard function Pi(q, omega) in the limit:")
    print(f"Momentum transfer q = {momentum_transfer_q}")
    print(f"External frequency omega = {external_frequency_omega}")

    print("\nIn this limit, the physical result is that the Lindhard function equals the")
    print("negative of the density of states at the Fermi energy, g(epsilon_F).")
    print("Pi(q=0, omega=0) = -g(epsilon_F)")

    print("\nTo obtain a universal, dimensionless numerical value, we normalize this result")
    print("by the density of states at the Fermi energy itself.")

    print("\nThe final equation for the normalized value is:")
    # Using the variable names to represent the physical quantities
    # The equation is: Pi(q, omega) / g(epsilon_F) = -1
    print(f"Pi(q={momentum_transfer_q}, omega={external_frequency_omega}) / g(epsilon_F) = {normalized_pi_value}")

    print("\n-------------------------------------------------")
    print(f"The requested numerical value is the result of this normalization.")

solve_lindhard_function()