import math

def analyze_beta_spectrum_setup():
    """
    Calculates the Larmor radius for the beta decay experiment parameters
    to verify the suitability of the proposed magnetic field.
    """
    # --- Constants ---
    Q = 1.60217663e-19  # Elementary charge in Coulombs
    M0_EV = 0.511e6       # Electron rest mass in eV
    C = 299792458       # Speed of light in m/s

    # --- Given Parameters ---
    kinetic_energy_MeV = 1.0
    magnetic_field_T = 166e-3 # 166 mT

    # Convert units for calculation
    kinetic_energy_J = kinetic_energy_MeV * 1e6 * Q
    m0c2_J = M0_EV * Q

    # --- Calculation ---
    # 1. Calculate total relativistic energy
    total_energy_J = kinetic_energy_J + m0c2_J

    # 2. Calculate relativistic momentum
    # E_total^2 = (p*c)^2 + (m0*c^2)^2  => p = sqrt(E_total^2 - (m0c2)^2) / c
    pc_J = math.sqrt(total_energy_J**2 - m0c2_J**2)
    momentum_kg_m_s = pc_J / C

    # 3. Calculate maximum Larmor (gyration) radius
    # r = p_perp / (q * B). Maximum radius occurs when perpendicular momentum p_perp is total momentum p.
    radius_m = momentum_kg_m_s / (Q * magnetic_field_T)

    # --- Output Results ---
    print("This script verifies that the proposed magnetic field is suitable for confining 1 MeV electrons.")
    print("\n--- Calculation of Maximum Larmor Radius ---")
    print(f"Kinetic Energy (K): {kinetic_energy_MeV} MeV")
    print(f"Magnetic Field (B): {magnetic_field_T * 1000} mT")
    print("-" * 40)

    print("Step 1: Calculate total momentum (p) of the electron.")
    print(f"   Relativistic momentum p = sqrt((K + m0*c^2)^2 - (m0*c^2)^2) / c")
    print(f"   p = {momentum_kg_m_s:.4e} kg*m/s")
    print()

    print("Step 2: Calculate the maximum Larmor radius (r_max).")
    print(f"   Larmor Radius r = p / (q * B)")
    # Final equation with numbers
    print(f"   r = {momentum_kg_m_s:.4e} / ({Q:.4e} * {magnetic_field_T})")
    print(f"   r = {radius_m:.4f} m")
    print(f"   r = {radius_m * 100:.2f} cm")
    print("-" * 40)

    print("\nConclusion:")
    print(f"The maximum gyration radius is {radius_m * 100:.2f} cm. This is a small radius, confirming that")
    print("a 166 mT field is strong enough to effectively confine the electrons and guide")
    print("them towards the detector in a typical laboratory setup.")

# Run the analysis
analyze_beta_spectrum_setup()