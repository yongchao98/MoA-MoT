import math

def calculate_gyroradius():
    """
    Calculates the maximum gyroradius for a beta particle in a magnetic field.

    This calculation demonstrates the feasibility of using a magnetic field to confine
    beta particles for spectroscopic measurements. A small gyroradius indicates
    that the particles will be tightly guided along the magnetic field lines.
    """
    # --- Constants ---
    # Kinetic energy of the electron in MeV
    KE_MeV = 1.0
    # Magnetic field strength in Tesla
    B_T = 166e-3
    # Electron rest mass in MeV/c^2
    m0_MeV_c2 = 0.511
    # Speed of light in m/s
    c_ms = 299792458
    # Elementary charge in Coulombs
    q_C = 1.602176634e-19
    # Conversion factor from MeV to Joules
    MeV_to_J = 1.602176634e-13

    # --- Calculation ---
    print("Step 1: Calculate the total relativistic energy (E) of the electron.")
    # Total energy E = Kinetic Energy (KE) + Rest Mass Energy (m0*c^2)
    E_MeV = KE_MeV + m0_MeV_c2
    print(f"E = KE + m0*c^2 = {KE_MeV:.1f} MeV + {m0_MeV_c2:.3f} MeV = {E_MeV:.3f} MeV")
    print("-" * 30)

    print("Step 2: Calculate the relativistic momentum (p) using the energy-momentum relation E^2 = (pc)^2 + (m0*c^2)^2.")
    # p = sqrt(E^2 - (m0*c^2)^2) / c
    p_MeV_c = math.sqrt(E_MeV**2 - m0_MeV_c2**2)
    print(f"p = sqrt(E^2 - (m0*c^2)^2) / c")
    print(f"p = sqrt({E_MeV:.3f}^2 - {m0_MeV_c2:.3f}^2) / c = {p_MeV_c:.3f} MeV/c")
    print("-" * 30)
    
    print("Step 3: Convert the momentum from MeV/c to SI units (kg*m/s).")
    # p [kg*m/s] = p [MeV/c] * (MeV_to_J / c)
    p_SI = p_MeV_c * (MeV_to_J / c_ms)
    print(f"p [kg*m/s] = {p_MeV_c:.3f} MeV/c * ({MeV_to_J:.6e} J/MeV / {c_ms:.0f} m/s) = {p_SI:.4e} kg*m/s")
    print("-" * 30)

    print("Step 4: Calculate the maximum gyroradius (r) using the formula r = p / (q * B).")
    # This assumes the momentum is entirely perpendicular to the field, giving the maximum possible radius.
    gyroradius_m = p_SI / (q_C * B_T)
    print("The final equation is: r = p / (q * B)")
    print(f"r = {p_SI:.4e} kg*m/s / ({q_C:.6e} C * {B_T:.3f} T)")
    
    gyroradius_cm = gyroradius_m * 100
    print(f"Final Result: The maximum gyroradius is {gyroradius_m:.4f} meters, or {gyroradius_cm:.2f} cm.")

if __name__ == '__main__':
    calculate_gyroradius()
