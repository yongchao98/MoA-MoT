import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # 1. Define physical constants and given parameters
    # High-precision constants in SI units
    hbar = 1.054571817e-34  # J·s
    amu_to_kg = 1.66053906660e-27 # kg
    e_charge = 1.602176634e-19 # C (for J to eV conversion)

    # Given parameters for H-Br
    r0 = 141.4e-12       # Bond length in meters
    k = 400.0             # Force constant in N/m
    mH_amu = 1.008        # Mass of Hydrogen in amu
    mBr_amu = 79.904      # Mass of Bromine in amu

    # 2. Calculate molecular properties
    # Convert masses to kg
    mH_kg = mH_amu * amu_to_kg
    mBr_kg = mBr_amu * amu_to_kg
    
    # Calculate reduced mass (μ) in kg
    mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)

    # Calculate moment of inertia (I) in kg·m^2
    I = mu * r0**2

    # 3. Calculate the centrifugal distortion constant (D_J) in Joules
    # Using the formula D_J = (ħ⁴ * μ) / (2 * I³ * k)
    D_Joule = (hbar**4 * mu) / (2 * I**3 * k)

    # 4. Calculate the energy shifts (ΔE) in Joules
    # For J=0 -> J=1 transition, shift is -4 * D_J
    shift_0_to_1_J = -4 * D_Joule
    
    # For J=1 -> J=2 transition, shift is -32 * D_J
    shift_1_to_2_J = -32 * D_Joule
    
    # 5. Convert energy shifts to quecto-electronvolts (qeV)
    # Conversion factor from Joules to qeV
    joule_to_qeV = 1.0 / (e_charge * 1e-30)

    shift_0_to_1_qeV = shift_0_to_1_J * joule_to_qeV
    shift_1_to_2_qeV = shift_1_to_2_J * joule_to_qeV

    # 6. Print the results in a detailed format
    print(f"Calculation of the centrifugal distortion constant D_J:")
    print(f"μ (reduced mass) = {mu:.4e} kg")
    print(f"I (moment of inertia) = {I:.4e} kg·m²")
    print(f"D_J = (ħ⁴ * μ) / (2 * I³ * k) = {D_Joule:.4e} J\n")

    print("1. Energy shift for the transition from J = 0 to J = 1:")
    print(f"   Formula: ΔE = -4 * D_J")
    print(f"   ΔE = -4 * {D_Joule:.4e} J")
    print(f"   ΔE = {shift_0_to_1_J:.4e} J")
    print(f"   ΔE = {shift_0_to_1_qeV:.4e} qeV\n")

    print("2. Energy shift for the transition from J = 1 to J = 2:")
    print(f"   Formula: ΔE = -4 * (1 + 1)³ * D_J = -32 * D_J")
    print(f"   ΔE = -32 * {D_Joule:.4e} J")
    print(f"   ΔE = {shift_1_to_2_J:.4e} J")
    print(f"   ΔE = {shift_1_to_2_qeV:.4e} qeV")

# Execute the function
calculate_energy_shifts()