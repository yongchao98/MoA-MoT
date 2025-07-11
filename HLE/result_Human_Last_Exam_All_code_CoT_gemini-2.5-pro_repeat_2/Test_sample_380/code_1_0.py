import math

def calculate_parker_spiral_angle():
    """
    Calculates the angle of the Parker Spiral at a given radial distance.

    This demonstrates that the interplanetary magnetic field at L1 (1 AU)
    is not purely radial, justifying why a field-aligned coordinate system
    is physically relevant, even if a radial system is used for convenience.
    """
    # Define constants
    # Solar angular velocity in radians per second
    omega_sun = 2.7e-6  # rad/s
    # Typical solar wind speed in meters per second (450 km/s)
    v_sw = 450 * 1000  # m/s
    # Radial distance for L1 point (1 Astronomical Unit) in meters
    r_au = 1.496e11  # m

    print("--- Parker Spiral Angle Calculation ---")
    print(f"This script calculates the angle between the interplanetary magnetic field and the radial direction at 1 AU.")
    print("\nInputs for the equation tan(ψ) = (Ω * r) / V_sw:")
    print(f"Solar Angular Velocity (Ω): {omega_sun:.2e} rad/s")
    print(f"Radial Distance (r): {r_au:.3e} m (1 AU)")
    print(f"Solar Wind Speed (V_sw): {v_sw:.0f} m/s")

    # Calculate the tangent of the angle
    tan_psi = (omega_sun * r_au) / v_sw
    
    # Calculate the angle in radians
    psi_rad = math.atan(tan_psi)
    
    # Convert the angle to degrees
    psi_deg = math.degrees(psi_rad)

    # The equation is tan(psi) = (omega_sun * r) / v_sw
    print("\nFinal Equation with values:")
    print(f"tan(ψ) = ({omega_sun:.2e} * {r_au:.3e}) / {v_sw:.0f}")
    print(f"tan(ψ) = {tan_psi:.4f}")
    print(f"ψ = atan({tan_psi:.4f})")
    print(f"ψ ≈ {psi_deg:.2f} degrees")
    
    print("\n--- Conclusion ---")
    print("As shown, the magnetic field at L1 is expected to be at an angle of ~42 degrees")
    print("to the radial direction. This confirms the field is not radial.")
    print("The use of a radial coordinate system for helicity is a practical choice for data analysis, not a physical assumption.")


if __name__ == '__main__':
    calculate_parker_spiral_angle()
<<<No, the local magnetic field at the L1 point is generally not radial due to the Parker Spiral effect. The use of a coordinate system with axes perpendicular to the radial direction is a standard convention for data analysis, prized for its simplicity and consistency, even though the physical wave propagation is along the local, non-radial magnetic field.>>>