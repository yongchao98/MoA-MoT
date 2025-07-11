import math

def explain_glucose_tube_appearance():
    """
    This script conceptually models why a tube of D-glucose solution with a light beam
    shone through it appears as a spiraling rainbow to a side observer.
    """

    print("--- Conceptual Analysis of Light in a D-Glucose Tube ---")
    print("1. Light Composition: White light is a mix of colors (wavelengths).")
    # We model three representative colors: Red, Green, and Blue.
    wavelengths = {"red": 650, "green": 550, "blue": 450} # in nanometers

    print("2. Optical Rotatory Dispersion (ORD):")
    print("D-glucose rotates the polarization of light. The angle of rotation depends on wavelength.")
    print("Shorter wavelengths (Blue) rotate MORE than longer wavelengths (Red).\n")
    
    # We'll use a simplified model for rotation angle: Angle = (k / wavelength^2) * distance
    # The 'k' value is arbitrary, chosen to show the effect clearly.
    # It incorporates concentration, specific rotation, etc.
    ROTATION_CONSTANT_K = 20000000  # An arbitrary constant for demonstration

    print("3. Observation from the Side:")
    print("An observer at the side sees scattered light. The color they see most clearly")
    print("depends on the polarization angle of the light at that point.\n")
    
    # Let's assume the observer sees a color best when its polarization has rotated by 90 degrees.
    # We will calculate the distance down the tube where this happens for each color.
    # Equation: 90 = (ROTATION_CONSTANT_K / wavelength^2) * distance
    # Solved for distance: distance = 90 * (wavelength^2) / ROTATION_CONSTANT_K
    
    print("--- Calculation of Viewing Distances (Observer at a fixed angle) ---")
    print("Distance at which each color's polarization rotates 90 degrees:")
    
    distances_for_90_deg = {}
    for color, wavelength_nm in wavelengths.items():
        # Equation solved for distance in meters
        distance_m = (90 * (wavelength_nm ** 2)) / ROTATION_CONSTANT_K
        distances_for_90_deg[color] = distance_m
        print(f"- {color.capitalize():<5}: For rotation = 90 deg, dist = {90} * ({wavelength_nm}^2) / {ROTATION_CONSTANT_K} = {distance_m:.3f} m")

    print("\nResult: The observer sees different colors at different distances, creating a rainbow along the tube.\n")
    
    print("--- Effect of Changing Viewing Angle ---")
    print("If the observer moves around the tube, the optimal angle for viewing changes.")
    print("Let's say the new optimal angle is 180 degrees.\n")
    
    # Equation: 180 = (ROTATION_CONSTANT_K / wavelength^2) * distance
    # Solved for distance: distance = 180 * (wavelength^2) / ROTATION_CONSTANT_K
    
    print("Distance at which each color's polarization rotates 180 degrees:")
    distances_for_180_deg = {}
    for color, wavelength_nm in wavelengths.items():
        # Equation solved for distance in meters
        distance_m = (180 * (wavelength_nm ** 2)) / ROTATION_CONSTANT_K
        distances_for_180_deg[color] = distance_m
        print(f"- {color.capitalize():<5}: For rotation = 180 deg, dist = {180} * ({wavelength_nm}^2) / {ROTATION_CONSTANT_K} = {distance_m:.3f} m")

    print("\nResult: The locations of the colors shift when the viewing angle changes.\n")
    
    print("--- Final Conclusion ---")
    print("The observed color depends on two variables:")
    print("  1. The distance along the tube (axial position).")
    print("  2. The viewing position around the tube (angular position).")
    print("A pattern that changes with both axial and angular position defines a SPIRAL or HELIX.")
    print("Therefore, the appearance is a rainbow gradient progressing in a spiral.")


# Run the explanation
explain_glucose_tube_appearance()