import numpy as np

def calculate_lindhard_function_limit():
    """
    Calculates the value of the Lindhard polarization function Pi(q, w)
    in the limit q->0 and w=0 for a 3D homogeneous electron gas.

    This script demonstrates that the dimensionless ratio Pi(0,0) / g(eps_F)
    is a universal constant, where g(eps_F) is the density of states at the
    Fermi level.
    """

    # Physical constants in SI units
    m_e = 9.1093837e-31      # Electron mass in kg
    hbar = 1.054571817e-34   # Reduced Planck constant in J*s

    # We choose an arbitrary electron density n (in m^-3) to show that the
    # final dimensionless ratio is independent of it.
    # Let's use a value typical for metals.
    n = 2.65e28

    # --- Calculation Steps ---

    # Step 1: Calculate the Fermi wavevector k_F from the electron density n.
    # The formula is k_F = (3 * pi^2 * n)^(1/3).
    k_F = (3 * np.pi**2 * n)**(1/3)

    # Step 2: Calculate the density of states at the Fermi level, g(eps_F).
    # This includes a factor of 2 for spin degeneracy.
    # The formula is g(eps_F) = (m_e * k_F) / (pi^2 * hbar^2).
    g_epsilon_F = (m_e * k_F) / (np.pi**2 * hbar**2)

    # Step 3: The Lindhard function at q=0, w=0 is the negative of g(eps_F).
    Pi_00 = -g_epsilon_F

    # Step 4: The universal dimensionless value is the ratio of Pi(0,0) to g(eps_F).
    dimensionless_value = Pi_00 / g_epsilon_F

    # --- Output ---
    print("The Lindhard function at zero momentum and frequency, Pi(0,0), is equal to the negative of the density of states at the Fermi level, g(eps_F).")
    print("The absolute value of Pi(0,0) depends on the electron density, but the ratio Pi(0,0) / g(eps_F) is a universal constant.")
    print("\nDemonstration with an example calculation:")
    print(f"For an electron density n = {n:.2e} m^-3:")
    print(f"  - The density of states g(eps_F) is calculated to be {g_epsilon_F:.3e} J^-1 m^-3.")
    print(f"  - Therefore, Pi(0,0) = -g(eps_F) = {Pi_00:.3e} J^-1 m^-3.")
    
    print("\nThe final equation for the universal dimensionless value is:")
    # The following line prints each number in the final equation as requested.
    print(f"Value = Pi(0,0) / g(eps_F) = ({Pi_00:.3e}) / ({g_epsilon_F:.3e}) = {dimensionless_value:.1f}")

if __name__ == '__main__':
    calculate_lindhard_function_limit()