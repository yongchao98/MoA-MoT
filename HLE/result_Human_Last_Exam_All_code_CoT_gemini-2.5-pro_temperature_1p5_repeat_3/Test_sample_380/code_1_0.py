import math

def calculate_parker_spiral_angle():
    """
    Calculates the Parker Spiral angle of the interplanetary magnetic field (IMF)
    at a given distance for a typical solar wind speed.

    The Parker Spiral angle (Psi) is the angle between the local magnetic field
    vector and the radial direction from the Sun.
    """

    # --- Constants ---
    # Sun's angular rotation speed in radians per second
    # (approx. 1 rotation every 27 days)
    omega_sun = 2 * math.pi / (27 * 24 * 3600)  # rad/s

    # Distance from Sun to L1 point (approx. 1 Astronomical Unit)
    r = 1.496e11  # meters

    # Typical solar wind speed in meters per second
    v_sw = 400 * 1000  # 400 km/s converted to m/s

    # --- Calculation ---
    # The formula is tan(Psi) = (omega_sun * r) / v_sw
    # Note: We take the absolute value as we are interested in the angle's magnitude.
    tan_psi = (omega_sun * r) / v_sw

    # Calculate the angle in radians using arctan
    psi_rad = math.atan(tan_psi)

    # Convert the angle to degrees for easier interpretation
    psi_deg = math.degrees(psi_rad)

    # --- Output ---
    print("--- Parker Spiral Angle Calculation at L1 (1 AU) ---")
    print(f"For a typical solar wind speed of {v_sw/1000:.0f} km/s:")
    print("The equation for the tangent of the spiral angle Psi is: tan(Psi) = (Omega * r) / V_sw")
    print(f"The calculated tangent value is: tan(Psi) = ({omega_sun:.2e} * {r:.3e}) / {v_sw:.1e} = {tan_psi:.4f}")
    print(f"This gives a Parker Spiral angle of {psi_deg:.2f} degrees from the radial direction.")
    print("\n--- Justification ---")
    print("No, the local magnetic field at the L1 point is not radial.")
    print("As the calculation shows, for a typical solar wind speed, the interplanetary magnetic field (IMF) makes an angle of approximately 45 degrees with the radial direction. This is known as the Parker Spiral.")
    print("\nThe reason magnetic helicity is often calculated using components perpendicular to the *radial* direction (e.g., Y and Z in GSE coordinates) is one of practical convenience and convention:")
    print("1. Standardized Coordinate System: Solar wind data is typically analyzed in a fixed coordinate system where the primary axis (X or R) points from the Sun to the spacecraft. This provides a stable, non-rotating frame of reference for statistical studies.")
    print("2. Dominant Flow Direction: The solar wind plasma flows almost perfectly radially away from the Sun. It is therefore natural to analyze wave fluctuations in the plane perpendicular to the main direction of flow and energy transport.")
    print("3. Approximation: While AIC waves propagate along the *local* magnetic field, using the radial direction as the reference axis is a simplification. It treats the dominant plasma flow direction as the primary axis, making the analysis simpler than constantly rotating into a field-aligned coordinate system, especially for long time-series data.")

if __name__ == '__main__':
    calculate_parker_spiral_angle()