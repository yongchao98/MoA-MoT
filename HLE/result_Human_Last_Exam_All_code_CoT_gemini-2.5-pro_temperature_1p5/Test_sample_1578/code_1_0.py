import math

def simulate_three_polarizers():
    """
    Simulates the Three-Polarizer Paradox to illustrate the principle in option C.

    The paradox:
    1. A vertical polarizer (P1 at 0 degrees) and a horizontal polarizer (P3 at 90 degrees)
       placed one after another will block all light.
    2. However, inserting a third polarizer (P2) between them at a 45-degree angle
       will unintuitively allow some light to pass through the entire system.

    This demonstrates that measurement (the middle polarizer) fundamentally changes the system,
    a principle that can be used to detect the demon's observation.
    """
    initial_intensity = 1.0  # Normalized initial unpolarized light intensity
    p1_angle = 0
    p2_angle = 45
    p3_angle = 90

    # According to Malus's Law, Intensity I = I_0 * cos^2(theta)
    # where theta is the angle between the light's polarization and the polarizer's axis.

    # After the first polarizer, intensity is halved and light is polarized at p1_angle
    intensity_after_p1 = initial_intensity / 2.0

    # Calculate the angle difference between P1 and P2
    theta_12 = math.radians(p2_angle - p1_angle)
    intensity_after_p2 = intensity_after_p1 * (math.cos(theta_12) ** 2)

    # Calculate the angle difference between P2 and P3
    theta_23 = math.radians(p3_angle - p2_angle)
    final_intensity = intensity_after_p2 * (math.cos(theta_23) ** 2)

    print("Simulating the Three-Polarizer Paradox:")
    print(f"Initial Light Intensity: {initial_intensity}")
    print(f"Polarizer 1 Angle: {p1_angle} degrees")
    print(f"Polarizer 2 Angle: {p2_angle} degrees")
    print(f"Polarizer 3 Angle: {p3_angle} degrees")
    print("-" * 30)

    # Without P2, the intensity would be 0
    intensity_without_p2 = intensity_after_p1 * (math.cos(math.radians(p3_angle - p1_angle)) ** 2)
    print(f"Result without middle polarizer (P2): {intensity_without_p2:.2f}")

    # With P2, some light passes
    print(f"Result with middle polarizer (P2): {final_intensity:.2f}")
    print("-" * 30)
    print("This happens because the middle polarizer 'measures' the photons and re-polarizes them,")
    print("allowing a component of the light to pass through the final polarizer.")
    print("\nFinal Equation Output:")
    # Per instructions, printing each number in the final equation:
    # Final Intensity = (Initial/2) * cos^2(theta_12) * cos^2(theta_23)
    # Final Intensity = 0.5 * cos^2(45) * cos^2(45) = 0.125
    # To simplify, we can treat the light after P1 as having intensity 1.
    # Then I_final = 1 * cos^2(45) * cos^2(45) = 0.5 * 0.5 = 0.25 (of what's left after P1)
    
    i1 = 1.0 # Let's say intensity after P1 is our new baseline
    angle1 = 45
    angle2 = 45 # (which is 90-45)
    i_final_simplified = i1 * (math.cos(math.radians(angle1))**2) * (math.cos(math.radians(angle2))**2)
    
    print(f"Based on a normalized intensity of {int(i1)} after the first polarizer:")
    print(f"Final Equation: I_final = {int(i1)} * cos^2({angle1}) * cos^2({angle2}) = {i_final_simplified:.3f}")


simulate_three_polarizers()