import scipy.constants as const

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define constants
    r0 = 141.4e-12      # Bond length in meters (m)
    k = 400.0           # Force constant in N/m
    mH_amu = 1.008      # Mass of Hydrogen in amu
    mBr_amu = 79.904    # Mass of Bromine in amu
    
    # Physical constants from scipy
    amu_to_kg = const.physical_constants['atomic mass unit-kilogram relationship'][0]
    hbar = const.hbar
    e = const.e

    print("--- Constants Used ---")
    print(f"Bond Length (r0): {r0} m")
    print(f"Force Constant (k): {k} N/m")
    print(f"Reduced Planck Constant (hbar): {hbar:.4e} J*s")
    print(f"Elementary Charge (e): {e:.4e} C\n")
    
    # 2. Calculate reduced mass (mu) in kg
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg
    mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)
    print("--- Intermediate Calculations ---")
    print(f"Reduced Mass (mu): {mu:.4e} kg")

    # 3. Calculate the centrifugal distortion constant (D_J) in Joules
    # Formula: D_J = hbar^4 / (2 * mu^2 * r0^6 * k)
    D_J = (hbar**4) / (2 * mu**2 * r0**6 * k)
    print(f"Centrifugal Distortion Constant (D_J): {D_J:.4e} J\n")

    # 4. Calculate energy shifts for each transition
    
    # For J=0 -> J=1
    print("--- Transition 1: J = 0 to J = 1 ---")
    factor1 = 4
    delta_E1_J = -factor1 * D_J
    print(f"Equation: Delta_E = -{factor1} * D_J")
    print(f"Calculation (Joules): Delta_E = -{factor1} * {D_J:.4e} J = {delta_E1_J:.4e} J")
    
    # For J=1 -> J=2
    print("\n--- Transition 2: J = 1 to J = 2 ---")
    factor2 = 32
    delta_E2_J = -factor2 * D_J
    print(f"Equation: Delta_E = -{factor2} * D_J")
    print(f"Calculation (Joules): Delta_E = -{factor2} * {D_J:.4e} J = {delta_E2_J:.4e} J")

    # 5. Convert from Joules to quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV = 1e-30 * e Joules
    joules_to_qeV = 1 / (e * 1e-30)
    delta_E1_qeV = delta_E1_J * joules_to_qeV
    delta_E2_qeV = delta_E2_J * joules_to_qeV

    print("\n--- Final Results in quecto-electronvolts (qeV) ---")
    print(f"Energy Shift for J=0 -> J=1: {delta_E1_qeV:.4e} qeV")
    print(f"Energy Shift for J=1 -> J=2: {delta_E2_qeV:.4e} qeV")

if __name__ == '__main__':
    calculate_energy_shifts()