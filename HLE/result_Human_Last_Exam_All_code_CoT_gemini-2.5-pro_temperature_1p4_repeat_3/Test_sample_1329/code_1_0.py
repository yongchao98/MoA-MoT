import math

def calculate_gyroradius():
    """
    Calculates the maximum gyroradius of a 1 MeV electron in a 166 mT magnetic field.
    This is to verify that the proposed field strength is physically reasonable for the experiment.
    """

    # --- Constants ---
    # Kinetic energy in Joules
    KE_eV = 1e6  # 1 MeV
    e_charge = 1.602e-19 # Elementary charge in Coulombs
    KE_J = KE_eV * e_charge

    # Electron rest mass in kg and MeV/c^2
    m0_kg = 9.109e-31
    m0_c2_eV = 0.511e6 # 0.511 MeV

    # Speed of light in m/s
    c = 2.998e8

    # Magnetic field in Tesla
    B = 0.166  # 166 mT

    print("--- Input Parameters ---")
    print(f"Kinetic Energy (KE): {KE_eV / 1e6} MeV")
    print(f"Magnetic Field (B): {B * 1000} mT")
    print("-" * 25)

    # --- Relativistic Calculation ---
    # Total energy E = KE + m0*c^2
    E_J = KE_J + (m0_c2_eV * e_charge)

    # From E^2 = (pc)^2 + (m0*c^2)^2, we find momentum p
    pc_squared = E_J**2 - (m0_c2_eV * e_charge)**2
    p = math.sqrt(pc_squared) / c

    print("--- Calculation Steps ---")
    print(f"Relativistic Momentum (p): {p:.4e} kg*m/s")

    # Gyroradius r = p / (q*B).
    # We use total momentum p for the worst-case scenario where velocity is fully perpendicular to B.
    q = e_charge
    r = p / (q * B)

    print(f"\nThe final equation for the maximum gyroradius (r) is: r = p / (q * B)")
    print(f"r = {p:.4e} kg*m/s / ({q:.4e} C * {B:.3f} T)")
    print(f"Calculated Maximum Gyroradius (r): {r:.4f} m or {r * 100:.2f} cm")
    print("-" * 25)

    print("\nConclusion:")
    print("A maximum gyroradius of a few centimeters is a manageable size for a laboratory vacuum chamber.")
    print("This confirms that 166 mT is a suitable field strength to guide and contain the 1 MeV beta particles.")

calculate_gyroradius()