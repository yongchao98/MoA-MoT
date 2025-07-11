import numpy as np

def calculate_parker_spiral_angle():
    """
    Calculates the Parker Spiral angle at 1 AU for a given solar wind speed.

    This demonstrates that the interplanetary magnetic field at the L1 point is
    not radial, justifying the need for careful coordinate system selection
    when calculating properties like magnetic helicity.
    """
    # Constants
    # Sun's angular velocity in radians per second
    omega_sun_rad_per_s = 2.9e-6
    # Distance from Sun (1 Astronomical Unit) in meters
    r_au_m = 1.5e11
    # Typical solar wind speed in meters per second (400 km/s)
    v_sw_m_per_s = 400 * 1000

    print("Calculating the Parker Spiral angle at the L1 point (~1 AU)...")
    print("-" * 60)
    print(f"Inputs:")
    print(f"  Sun's Angular Velocity (Ω): {omega_sun_rad_per_s:.2e} rad/s")
    print(f"  Radial Distance (R): {r_au_m:.2e} m")
    print(f"  Solar Wind Speed (V_sw): {v_sw_m_per_s:.1f} m/s (or {v_sw_m_per_s/1000:.0f} km/s)")
    print("-" * 60)

    # Calculate tan(psi)
    tan_psi = (omega_sun_rad_per_s * r_au_m) / v_sw_m_per_s
    # Calculate the angle in radians
    psi_rad = np.arctan(tan_psi)
    # Convert angle to degrees
    psi_deg = np.degrees(psi_rad)

    # Print the equation with values
    print("Equation: tan(ψ) = (Ω * R) / V_sw")
    print(f"Calculation: tan(ψ) = ({omega_sun_rad_per_s:.2e} * {r_au_m:.2e}) / {v_sw_m_per_s:.1f}")
    print(f"Resulting tan(ψ) = {tan_psi:.4f}")
    print("-" * 60)

    print(f"The calculated Parker Spiral angle (ψ) is: {psi_deg:.2f} degrees.")
    print("\nConclusion:")
    print("Since this angle is significantly different from 0 degrees, the local magnetic field")
    print("at the L1 point is NOT radial. Therefore, using magnetic field components")
    print("perpendicular to the radial direction to calculate helicity is an approximation.")


calculate_parker_spiral_angle()
<<<The use of magnetic field components perpendicular to the radial direction is a simplifying approximation. The local magnetic field at the L1 point is typically not radial due to the Parker Spiral effect, making an angle of about 45 degrees with the radial direction for average solar wind speeds. The correct but more complex method is to use a field-aligned coordinate system, often found using Minimum Variance Analysis.>>>