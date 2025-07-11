import math
import scipy.constants as const

def calculate_centrifugal_distortion_shifts():
    """
    Calculates the energy shifts due to centrifugal distortion for H-Br
    for J=0->1 and J=1->2 transitions.
    """
    # Step 1: Define given constants and fundamental constants
    # Given values
    r0_pm = 141.4    # bond length in picometers
    k_Nm = 400.0     # force constant in N/m
    mH_amu = 1.008   # mass of Hydrogen in amu
    mBr_amu = 79.904 # mass of Bromine in amu

    # Convert to SI units
    r0 = r0_pm * 1e-12    # m
    mH_kg = mH_amu * const.u # kg
    mBr_kg = mBr_amu * const.u # kg

    # Fundamental constants
    h = const.h
    c = const.c
    e = const.e
    pi = const.pi

    # Step 2: Calculate intermediate molecular properties
    # Reduced mass (mu) in kg
    mu = (mH_kg * mBr_kg) / (mH_kg + mBr_kg)
    # Moment of inertia (I) in kg*m^2
    I = mu * r0**2

    # Step 3: Calculate spectroscopic constants in SI wavenumbers (m^-1)
    # Rotational constant B in m^-1
    B_wavenumber = h / (8 * pi**2 * c * I)
    # Vibrational frequency omega_e in m^-1
    omega_e_wavenumber = (1 / (2 * pi * c)) * math.sqrt(k_Nm / mu)
    # Centrifugal distortion constant D in m^-1, from D = 4B^3 / omega_e^2
    D_wavenumber = (4 * B_wavenumber**3) / (omega_e_wavenumber**2)

    # Step 4 & 5: Calculate energy shifts and convert to qeV
    # The energy shift (ΔE) for a transition J_i -> J_f is:
    # ΔE = h*c*D * [ J_i^2*(J_i+1)^2 - J_f^2*(J_f+1)^2 ]
    # The final answer is requested in quecto-electronvolts (qeV)
    # 1 qeV = 1e-30 eV = 1e-30 * e Joules
    joules_to_qeV = 1.0 / (e * 1e-30)

    # --- Transition 1: J=0 to J=1 ---
    print("For the rotational transition from J = 0 to J = 1:")
    # The specific formula for this shift is: ΔE = h*c*D * [0^2(1)^2 - 1^2(2)^2]
    # This simplifies to ΔE = -4 * h * c * D
    factor1 = -4.0
    delta_E_J1 = factor1 * h * c * D_wavenumber
    delta_E_qeV1 = delta_E_J1 * joules_to_qeV
    
    print(f"The energy shift equation is: ΔE = -4 * h * c * D")
    print(f"ΔE = -4 * ({h:.4e} J·s) * ({c:.4e} m/s) * ({D_wavenumber:.4e} m⁻¹)")
    print(f"The final energy shift is {delta_E_qeV1:.4f} qeV")
    print("-" * 50)

    # --- Transition 2: J=1 to J=2 ---
    print("For the rotational transition from J = 1 to J = 2:")
    # The specific formula is: ΔE = h*c*D * [1^2(2)^2 - 2^2(3)^2]
    # This simplifies to ΔE = -32 * h * c * D
    factor2 = -32.0
    delta_E_J2 = factor2 * h * c * D_wavenumber
    delta_E_qeV2 = delta_E_J2 * joules_to_qeV

    print(f"The energy shift equation is: ΔE = -32 * h * c * D")
    print(f"ΔE = -32 * ({h:.4e} J·s) * ({c:.4e} m/s) * ({D_wavenumber:.4e} m⁻¹)")
    print(f"The final energy shift is {delta_E_qeV2:.4f} qeV")

if __name__ == "__main__":
    calculate_centrifugal_distortion_shifts()