import math

def analyze_moon_orbit():
    """
    Calculates the ratio of the Sun's gravitational force on the Moon to the
    Earth's gravitational force on the Moon to determine the shape of the
    Moon's orbit in a heliocentric frame.
    """
    # Physical constants
    M_sun = 1.989e30  # Mass of the Sun in kg
    M_earth = 5.972e24 # Mass of the Earth in kg
    # The mass of the moon is not needed as it cancels out in the ratio.
    
    # Average distances in meters
    # The distance from the Sun to the Moon is approximately the same as the Sun to the Earth
    R_sun_moon = 1.496e11 
    R_earth_moon = 3.844e8
    
    # The gravitational force formula is F = G * M1 * M2 / R^2.
    # When calculating the ratio of forces on the Moon, G and the Moon's mass cancel out.
    # Ratio = (F_sun_on_moon) / (F_earth_on_moon)
    # Ratio = (M_sun / R_sun_moon^2) / (M_earth / R_earth_moon^2)
    
    force_component_sun = M_sun / (R_sun_moon**2)
    force_component_earth = M_earth / (R_earth_moon**2)
    
    ratio = force_component_sun / force_component_earth
    
    print("To understand the Moon's path around the Sun, let's compare the gravitational forces acting on it.")
    print(f"The Sun's mass is {M_sun:.3e} kg.")
    print(f"The Earth's mass is {M_earth:.3e} kg.")
    print(f"The Sun-Moon distance is approximately {R_sun_moon:.3e} m.")
    print(f"The Earth-Moon distance is {R_earth_moon:.3e} m.")
    print("\nThe ratio of the Sun's gravitational pull to the Earth's gravitational pull on the Moon is:")
    print(f"(Mass_Sun / Dist_Sun_Moon^2) / (Mass_Earth / Dist_Earth_Moon^2)")
    print(f"= ({M_sun:.3e} / {R_sun_moon:.3e}^2) / ({M_earth:.3e} / {R_earth_moon:.3e}^2)")
    print(f"= {ratio:.2f}")

    print("\nConclusion:")
    print(f"The Sun's gravitational force on the Moon is more than {math.floor(ratio)} times stronger than the Earth's force.")
    print("This means the Moon's trajectory is always being pulled towards the Sun.")
    print("Therefore, the Moon's path is always concave towards the Sun, like a wave that never curves back on itself.")
    print("Looking at the options:")
    print(" - A is a straight line, incorrect.")
    print(" - B is a smooth curve without the Moon's wobble, incorrect.")
    print(" - D and E show loops/cusps where the path becomes convex (curves away from the Sun), which is physically impossible.")
    print(" - C correctly shows a wavy path that is always concave, representing the Moon orbiting the Earth while both are pulled strongly into orbit around the Sun.")
    print("\nThe correct choice is C.")

analyze_moon_orbit()