import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # --- 1. Define Constants ---
    # Given values
    m_H_amu = 1.008       # amu
    m_Br_amu = 79.904      # amu
    r0_pm = 141.4         # pm
    k_Nm = 400.0          # N/m

    # Fundamental constants
    u = 1.66053906660e-27 # kg/amu (atomic mass unit)
    h = 6.62607015e-34    # J*s (Planck constant)
    c_ms = 299792458.0    # m/s (speed of light in m/s)
    e = 1.602176634e-19   # C (elementary charge for J to eV conversion)
    
    # --- 2. Calculate Molecular Properties in SI units ---
    # Convert inputs to SI
    r0_m = r0_pm * 1e-12
    m_H_kg = m_H_amu * u
    m_Br_kg = m_Br_amu * u
    c_cms = c_ms * 100 # speed of light in cm/s for wavenumber calculations
    
    # Calculate reduced mass (mu)
    mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

    # Calculate moment of inertia (I)
    I = mu_kg * r0_m**2
    
    # --- 3. Calculate Spectroscopic Constants (in cm^-1) ---
    # Rotational constant B_bar in cm^-1
    B_bar_cm = h / (8 * math.pi**2 * c_cms * I)
    
    # Vibrational frequency omega_e_bar in cm^-1
    omega_e_bar_cm = (1 / (2 * math.pi * c_cms)) * math.sqrt(k_Nm / mu_kg)
    
    # Centrifugal distortion constant D_bar in cm^-1
    D_bar_cm = (4 * B_bar_cm**3) / (omega_e_bar_cm**2)

    # --- 4. Calculate Energy Shifts ---
    # Convert D from cm^-1 to Joules
    # Energy (J) = wavenumber (cm^-1) * 100 * h * c_ms
    D_J = D_bar_cm * 100 * h * c_ms

    print(f"The centrifugal distortion constant D is {D_J:.4e} Joules.")
    print("-" * 50)
    
    # Transition 1: J = 0 to J = 1
    J1 = 0
    shift1_J = -4 * D_J * (J1 + 1)**3
    shift1_eV = shift1_J / e
    shift1_qeV = shift1_eV * 1e30

    print("1. For the rotational transition from J = 0 to J = 1:")
    print(f"   ΔE = -4 * D * (J+1)³")
    print(f"   ΔE = -4 * {D_J:.4e} J * ({J1}+1)³ = {shift1_J:.4e} J")
    print(f"   In quecto-electronvolts, the energy shift is {shift1_qeV:.4f} qeV.")
    print("-" * 50)

    # Transition 2: J = 1 to J = 2
    J2 = 1
    shift2_J = -4 * D_J * (J2 + 1)**3
    shift2_eV = shift2_J / e
    shift2_qeV = shift2_eV * 1e30

    print("2. For the rotational transition from J = 1 to J = 2:")
    print(f"   ΔE = -4 * D * (J+1)³")
    print(f"   ΔE = -4 * {D_J:.4e} J * ({J2}+1)³ = {shift2_J:.4e} J")
    print(f"   In quecto-electronvolts, the energy shift is {shift2_qeV:.4f} qeV.")

calculate_energy_shifts()