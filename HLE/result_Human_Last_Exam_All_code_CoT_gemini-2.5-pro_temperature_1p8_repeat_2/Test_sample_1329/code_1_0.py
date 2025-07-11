import math
import scipy.constants as const

def calculate_gyroradius():
    """
    Calculates the gyroradius for a 1 MeV electron in a 166 mT magnetic field.
    This serves to verify the plausibility of the physical parameters given in the problem.
    """

    # --- Given Parameters ---
    kinetic_energy_MeV = 1.0  # Maximum kinetic energy in Mega-electron Volts
    magnetic_field_T = 166e-3  # Magnetic field in Tesla (166 mT)

    # --- Physical Constants ---
    # Electron rest energy in MeV
    electron_rest_energy_MeV = (const.m_e * const.c**2) / const.e / 1e6
    # Electron charge in Coulombs
    electron_charge_C = const.e

    # --- Calculations ---
    # 1. Total energy E = T + m_e*c^2
    total_energy_MeV = kinetic_energy_MeV + electron_rest_energy_MeV

    # 2. Relativistic momentum from E^2 = (pc)^2 + (m_e*c^2)^2
    #    pc = sqrt(E^2 - (m_e*c^2)^2)
    pc_MeV = math.sqrt(total_energy_MeV**2 - electron_rest_energy_MeV**2)

    # 3. Convert momentum from MeV/c to SI units (kg*m/s)
    #    p = (pc_MeV * 10^6 * e) / c
    momentum_SI = (pc_MeV * 1e6 * const.e) / const.c

    # 4. Gyroradius r = p / (q*B)
    gyroradius_m = momentum_SI / (electron_charge_C * magnetic_field_T)

    # --- Print Results ---
    print("Calculation to verify the physical scale of the problem:")
    print("-" * 50)
    print(f"Given maximum kinetic energy (T): {kinetic_energy_MeV:.3f} MeV")
    print(f"Given magnetic field (B): {magnetic_field_T * 1000:.0f} mT")
    print(f"Electron rest energy (m_e*c^2): {electron_rest_energy_MeV:.3f} MeV")
    print("-" * 50)
    
    print("Final Equation: r = p / (q * B)")
    print(f"The calculated gyroradius for the most energetic electrons is:")
    # The instructions require printing the numbers in the final equation.
    # r = p / (q * B)
    print(f"r = {momentum_SI:.3e} kg*m/s / ({electron_charge_C:.3e} C * {magnetic_field_T:.3f} T)")
    print(f"r = {gyroradius_m:.4f} m or {gyroradius_m * 100:.2f} cm")
    print("-" * 50)
    print("This radius is a few centimeters, which is a reasonable size for confining")
    print("electrons within a typical laboratory vacuum chamber.")
    print("Therefore, the given field strength is appropriate for the energy scale.")


# Run the calculation
calculate_gyroradius()