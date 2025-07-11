import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Define Physical and Given Constants ---

    # Given parameters for H-Br molecule
    r0_pm = 141.4  # Bond length in picometers
    k_Nm = 400.0   # Force constant in N/m
    m_H_amu = 1.008  # Mass of Hydrogen in amu
    m_Br_amu = 79.904 # Mass of Bromine in amu

    # Physical constants
    hbar = 1.054571817e-34  # Reduced Planck constant in J·s
    amu_to_kg = 1.66053906660e-27 # Conversion factor from amu to kg
    pm_to_m = 1e-12             # Conversion factor from pm to m
    e_charge = 1.602176634e-19    # Elementary charge in Coulombs (for J to eV conversion)
    eV_to_qeV = 1e-30           # Conversion factor from eV to qeV

    # --- 2. Convert to SI Units ---

    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    r0_m = r0_pm * pm_to_m

    # --- 3. Calculate Reduced Mass (mu) ---

    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # --- 4. Calculate Centrifugal Distortion Constant (D) in Joules ---
    
    # Using the formula D = hbar^4 / (2 * mu^2 * r0^6 * k)
    D_joules = (hbar**4) / (2 * mu_kg**2 * r0_m**6 * k_Nm)

    # --- 5. Calculate Energy Shifts (Delta_E) in Joules ---
    
    # For transition J=0 to J=1, the initial state is J=0
    J1 = 0
    delta_E1_joules = -4 * D_joules * (J1 + 1)**3
    
    # For transition J=1 to J=2, the initial state is J=1
    J2 = 1
    delta_E2_joules = -4 * D_joules * (J2 + 1)**3

    # --- 6. Convert Energy Shifts to quecto-electronvolts (qeV) ---
    
    # Conversion factor from Joules to qeV
    joules_to_qeV = 1 / (e_charge * eV_to_qeV)
    
    delta_E1_qeV = delta_E1_joules * joules_to_qeV
    delta_E2_qeV = delta_E2_joules * joules_to_qeV

    # --- 7. Print the Results ---
    print("Calculation for H-Br molecule:")
    print("-" * 30)
    print(f"Centrifugal distortion constant D = {D_joules:.4e} J")
    print("-" * 30)

    print("Energy shift for the J=0 to J=1 transition:")
    print(f"ΔE(J=0→1) = -4 * D * (0+1)^3")
    print(f"ΔE(J=0→1) = {delta_E1_qeV:.4e} qeV")
    print("")

    print("Energy shift for the J=1 to J=2 transition:")
    print(f"ΔE(J=1→2) = -4 * D * (1+1)^3")
    print(f"ΔE(J=1→2) = {delta_E2_qeV:.4e} qeV")

calculate_energy_shifts()