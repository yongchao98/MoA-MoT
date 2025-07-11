import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Define Constants ---
    # Given molecular properties
    r0_pm = 141.4  # Bond length in picometers
    k_Nm = 400.0   # Force constant in N/m
    mH_amu = 1.008  # Mass of Hydrogen in amu
    mBr_amu = 79.904 # Mass of Bromine in amu

    # Physical constants
    h = 6.62607015e-34      # Planck constant in J·s
    hbar = h / (2 * math.pi) # Reduced Planck constant in J·s
    amu_to_kg = 1.660539e-27 # Conversion factor from amu to kg
    pm_to_m = 1e-12          # Conversion factor from pm to m
    J_per_eV = 1.602176634e-19 # Joules per electronvolt
    qeV_per_eV = 1e30        # quecto-electronvolts per electronvolt

    # --- 2. Calculate Molecular Parameters in SI Units ---
    # Convert given properties to SI units
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg
    r0_m = r0_pm * pm_to_m

    # Calculate reduced mass (mu)
    mu_kg = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)

    # Calculate moment of inertia (I)
    I = mu_kg * r0_m**2

    # Calculate rotational constant (B) in Joules
    B_J = (hbar**2) / (2 * I)

    # Calculate angular vibrational frequency (omega) in rad/s
    omega_rad_s = math.sqrt(k_Nm / mu_kg)
    
    # Calculate the fundamental vibrational energy (omega_e) in Joules
    omega_e_J = hbar * omega_rad_s

    # Calculate centrifugal distortion constant (D) in Joules
    D_J = (4 * B_J**3) / (omega_e_J**2)
    
    # Conversion factor from Joules to quecto-electronvolts
    J_to_qeV = (1 / J_per_eV) * qeV_per_eV

    # --- 3. Calculate and Print Energy Shifts ---

    # For the transition from J = 0 to J = 1
    J1 = 0
    shift1_J = -4 * D_J * (J1 + 1)**3
    shift1_qeV = shift1_J * J_to_qeV

    print(f"For the transition from J = {J1} to J = {J1 + 1}:")
    print(f"The energy shift is calculated using the formula: ΔE = -4 * D * (J_initial + 1)³")
    print(f"   D = {D_J:.5e} J")
    print(f"   J_initial = {J1}")
    print(f"Equation: ΔE = -4 * ({D_J:.5e}) * ({J1} + 1)³ = {shift1_J:.5e} J")
    print(f"The energy shift is {shift1_qeV:.4e} qeV.\n")

    # For the transition from J = 1 to J = 2
    J2 = 1
    shift2_J = -4 * D_J * (J2 + 1)**3
    shift2_qeV = shift2_J * J_to_qeV

    print(f"For the transition from J = {J2} to J = {J2 + 1}:")
    print(f"The energy shift is calculated using the formula: ΔE = -4 * D * (J_initial + 1)³")
    print(f"   D = {D_J:.5e} J")
    print(f"   J_initial = {J2}")
    print(f"Equation: ΔE = -4 * ({D_J:.5e}) * ({J2} + 1)³ = {shift2_J:.5e} J")
    print(f"The energy shift is {shift2_qeV:.4e} qeV.")


# Execute the function to print the results
calculate_energy_shifts()