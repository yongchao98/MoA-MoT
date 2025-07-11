import math

def calculate_energy_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br rotational transitions.
    """
    # Step 1: Define constants and given parameters in SI units
    m_H_amu = 1.008      # mass of hydrogen in amu
    m_Br_amu = 79.904    # mass of bromine in amu
    r0 = 141.4e-12       # bond length in meters (pm -> m)
    k = 400              # force constant in N/m
    
    # Physical constants
    h = 6.62607015e-34   # Planck's constant in J*s
    amu_to_kg = 1.660539e-27 # conversion factor for amu to kg
    eV_to_J = 1.602176634e-19 # conversion for eV to Joules
    
    # Step 2: Calculate reduced mass (mu) in kg
    m_H_kg = m_H_amu * amu_to_kg
    m_Br_kg = m_Br_amu * amu_to_kg
    mu = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)
    
    # Step 3: Calculate the rotational constant B and distortion constant D
    # Moment of inertia (I)
    I = mu * r0**2
    
    # Rotational constant B in Hz
    B_hz = h / (8 * math.pi**2 * I)
    
    # Vibrational frequency (nu) in Hz
    nu_hz = (1 / (2 * math.pi)) * math.sqrt(k / mu)
    
    # Centrifugal distortion constant D in Hz
    D_hz = (4 * B_hz**3) / (nu_hz**2)
    
    # Calculate the energy term hD in Joules
    hD_joules = h * D_hz
    
    # Step 4: Calculate energy shift for each transition
    
    # Transition 1: J = 0 to J = 1
    J_initial_1 = 0
    # Final equation for the energy shift is: Delta_E = -4 * hD * (J_initial + 1)^3
    shift_joules_1 = -4 * hD_joules * (J_initial_1 + 1)**3
    
    # Convert Joules to qeV (1 qeV = 1e-30 eV)
    qeV_conversion_factor = eV_to_J * 1e-30
    shift_qeV_1 = shift_joules_1 / qeV_conversion_factor
    
    print("For the J=0 to J=1 transition:")
    print(f"The final equation for the shift is: -4 * (h*D) * ({J_initial_1} + 1)^3")
    print(f"Calculated shift = {shift_qeV_1:.4e} qeV")
    print("-" * 30)

    # Transition 2: J = 1 to J = 2
    J_initial_2 = 1
    # Final equation for the energy shift is: Delta_E = -4 * hD * (J_initial + 1)^3
    shift_joules_2 = -4 * hD_joules * (J_initial_2 + 1)**3
    
    # Convert Joules to qeV
    shift_qeV_2 = shift_joules_2 / qeV_conversion_factor
    
    print("For the J=1 to J=2 transition:")
    print(f"The final equation for the shift is: -4 * (h*D) * ({J_initial_2} + 1)^3")
    print(f"Calculated shift = {shift_qeV_2:.4e} qeV")


if __name__ == '__main__':
    calculate_energy_shifts()