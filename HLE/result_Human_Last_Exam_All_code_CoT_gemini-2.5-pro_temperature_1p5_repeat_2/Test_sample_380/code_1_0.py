import math

def calculate_parker_spiral_angle():
    """
    Calculates the Parker Spiral angle of the interplanetary magnetic field (IMF)
    at 1 AU for a typical solar wind speed.
    """
    # --- Constants ---
    # Sun's angular rotation speed in radians per second
    # (approx. for a 27-day sidereal rotation period)
    omega_sun = 2 * math.pi / (27 * 24 * 3600)  # rad/s

    # Distance from Sun (1 Astronomical Unit) in meters
    r_au = 1.496e11  # m

    # Typical solar wind speed in meters per second
    v_sw = 400e3  # m/s (400 km/s)

    print("--- Parker Spiral Angle Calculation ---")
    print(f"Using the following values for L1 (approx. 1 AU):")
    print(f"Sun's angular velocity (Ω): {omega_sun:.2e} rad/s")
    print(f"Radial distance (r): {r_au:.3e} m")
    print(f"Solar wind speed (V_sw): {v_sw:.1e} m/s")
    print("-" * 35)

    # --- Calculation ---
    # The equation for the angle Psi is: tan(Psi) = (Ω * r) / V_sw
    # Calculate the tangential speed due to Sun's rotation at distance r
    v_rotation_at_r = omega_sun * r_au
    
    # Calculate the ratio for tan(Psi)
    tan_psi = v_rotation_at_r / v_sw

    # Calculate the angle in radians and then convert to degrees
    psi_rad = math.atan(tan_psi)
    psi_deg = math.degrees(psi_rad)

    # --- Output the Final "Equation" and Result ---
    # The prompt asks to output each number in the final equation.
    print("Final Equation Steps:")
    print(f"Tangential speed at r (Ω * r) = {omega_sun:.2e} rad/s * {r_au:.3e} m = {v_rotation_at_r:.3f} m/s")
    print(f"tan(Ψ) = (Tangential speed) / (Solar wind speed) = {v_rotation_at_r:.3f} m/s / {v_sw:.1e} m/s = {tan_psi:.4f}")
    print(f"Ψ = atan({tan_psi:.4f})")
    print(f"\nThe calculated Parker Spiral angle (Ψ) is: {psi_deg:.2f} degrees")
    print("-" * 35)
    print("\nThis calculation shows the average magnetic field is not radial, but at an angle of ~45 degrees to the Sun-Earth line.")


calculate_parker_spiral_angle()
