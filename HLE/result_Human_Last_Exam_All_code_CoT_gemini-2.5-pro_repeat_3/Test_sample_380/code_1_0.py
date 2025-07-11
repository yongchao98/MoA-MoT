import numpy as np

def calculate_parker_spiral_angle():
    """
    Calculates the angle of the interplanetary magnetic field (Parker Spiral)
    with respect to the radial direction at 1 AU.
    """
    # Define constants
    # Sun's sidereal rotation period is ~25.38 days
    T_sun = 25.38 * 24 * 3600  # Sun's rotation period in seconds
    omega_sun = 2 * np.pi / T_sun  # Sun's angular velocity in rad/s

    # Distance from Sun (1 Astronomical Unit)
    r = 1.496e11  # in meters

    # Typical solar wind speed
    V_sw = 400 * 1000  # 400 km/s, converted to m/s

    # The formula for the Parker spiral angle (psi) is tan(psi) = (omega * r) / V_sw
    # We calculate the argument of arctan first
    tan_psi_arg = (omega_sun * r) / V_sw
    
    # Calculate the angle in radians
    angle_rad = np.arctan(tan_psi_arg)

    # Convert the angle to degrees for better intuition
    angle_deg = np.rad2deg(angle_rad)

    print("--- Parker Spiral Angle Calculation at 1 AU ---")
    print(f"This calculation shows why the magnetic field at L1 is not radial.")
    print("\nInputs for the equation: angle = arctan((Ω * r) / V_sw)")
    print(f"  Sun's angular velocity (Ω): {omega_sun:.2e} rad/s")
    print(f"  Radial distance (r): {r:.3e} m (1 AU)")
    print(f"  Solar wind speed (V_sw): {V_sw / 1000:.0f} km/s")

    print(f"\nCalculation steps:")
    print(f"  Numerator (Ω * r): {omega_sun * r:.3e} m/s")
    print(f"  Denominator (V_sw): {V_sw:.3e} m/s")
    print(f"  Argument for arctan ((Ω * r) / V_sw): {tan_psi_arg:.4f}")
    
    print(f"\nFinal Result:")
    print(f"The Parker spiral angle at 1 AU is {angle_deg:.2f} degrees from the radial direction.")
    print("\nConclusion: Since the angle is significantly different from 0, the local magnetic field is not radial.")
    print("Therefore, a coordinate transformation is necessary to use components perpendicular to the local field for helicity calculations.")

if __name__ == '__main__':
    calculate_parker_spiral_angle()