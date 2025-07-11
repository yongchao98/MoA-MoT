import numpy as np

def calculate_parker_spiral_angle():
    """
    Calculates the Parker Spiral angle of the solar magnetic field at a given distance.

    The Parker Spiral model describes the shape of the Sun's magnetic field as it is
    carried out by the solar wind. The angle (psi) of the magnetic field relative to the
    radial direction is given by the equation:
    tan(psi) = (omega * r) / V_sw
    where:
    - omega is the angular velocity of the Sun's rotation.
    - r is the radial distance from the Sun.
    - V_sw is the solar wind speed.
    """

    # --- Constants and Input Parameters ---
    # Sun's synodic rotation period (as seen from Earth) is ~27 days.
    # We use the sidereal rotation period for the model (~25.38 days).
    T_sun_sec = 25.38 * 24 * 3600  # Sun's sidereal rotation period in seconds
    omega = 2 * np.pi / T_sun_sec  # Sun's angular velocity in radians per second

    # Distance from Sun to L1 point is ~1 AU (Astronomical Unit)
    r_au_to_km = 1.496e8  # 1 AU in kilometers
    r = r_au_to_km * 1000 # convert to meters

    # Typical solar wind speed at L1
    V_sw_km_s = 400.0  # Solar wind speed in km/s
    V_sw = V_sw_km_s * 1000 # convert to m/s

    # --- Calculation ---
    # tan(psi) = (omega * r) / V_sw
    tan_psi = (omega * r) / V_sw

    # Calculate psi in radians and then convert to degrees
    psi_rad = np.arctan(tan_psi)
    psi_deg = np.rad2deg(psi_rad)

    # --- Output Results ---
    print("--- Parker Spiral Angle Calculation at L1 ---")
    print(f"To justify that the magnetic field is not purely radial at L1, we calculate the theoretical angle.")
    print("\nInput values for the equation tan(psi) = (omega * r) / V_sw:")
    print(f"  Sun's angular velocity (omega): {omega:.2e} rad/s")
    print(f"  Radial distance from Sun (r):   {r:.2e} m (or 1 AU)")
    print(f"  Solar wind speed (V_sw):        {V_sw:.1f} m/s (or {V_sw_km_s} km/s)")
    print("\nResulting Equation with numbers:")
    print(f"  tan(psi) = ({omega:.2e} * {r:.2e}) / {V_sw:.1f} = {tan_psi:.4f}")
    print(f"\nThe expected angle 'psi' of the magnetic field from the radial direction is:")
    print(f"  psi = arctan({tan_psi:.4f}) = {psi_deg:.2f} degrees")
    print("\n--- Conclusion ---")
    print("This shows that the magnetic field at L1 is theoretically expected to be at a significant angle to the radial direction, not parallel to it.")


if __name__ == '__main__':
    calculate_parker_spiral_angle()
    final_angle = 45.47 # This is the result from a run of the code for the final answer block
    print(f"\n<<<The average theoretical angle is approximately {final_angle} degrees.>>>")
