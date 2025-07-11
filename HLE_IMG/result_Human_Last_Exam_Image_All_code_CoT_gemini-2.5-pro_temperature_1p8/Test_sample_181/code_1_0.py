import math

def analyze_moon_orbit():
    """
    Calculates the ratio of the Sun's gravitational pull to the Earth's
    gravitational pull on the Moon to determine the shape of its orbit in an
    inertial frame centered on the Sun.
    """
    # Physical constants
    M_sun = 1.989e30  # Mass of the Sun in kg
    M_earth = 5.972e24 # Mass of the Earth in kg
    d_sun_moon = 1.496e11 # Average distance from the Sun to the Moon (approx. 1 AU) in meters
    d_earth_moon = 3.844e8  # Average distance from the Earth to the Moon in meters

    # Gravitational force is proportional to M / d^2
    # The gravitational constant G and the Moon's mass cancel out in the ratio.
    
    # Proportional force from the Sun on the Moon
    force_sun_prop = M_sun / (d_sun_moon**2)
    
    # Proportional force from the Earth on the Moon
    force_earth_prop = M_earth / (d_earth_moon**2)
    
    # Ratio of the two forces
    force_ratio = force_sun_prop / force_earth_prop
    
    print("This script determines the shape of the Moon's orbit around the Sun.")
    print("It does this by comparing the gravitational force from the Sun to the force from the Earth.")
    print("\nConstants used:")
    print(f"Mass of Sun: {M_sun:.3e} kg")
    print(f"Mass of Earth: {M_earth:.3e} kg")
    print(f"Sun-Moon distance: {d_sun_moon:.3e} m")
    print(f"Earth-Moon distance: {d_earth_moon:.3e} m")

    print(f"\nCalculation: Ratio = (M_sun / d_sun_moon^2) / (M_earth / d_earth_moon^2)")
    print(f"Ratio = ({M_sun:.3e} / {d_sun_moon:.3e}^2) / ({M_earth:.3e} / {d_earth_moon:.3e}^2)")
    print(f"Ratio = {force_sun_prop:.3e} / {force_earth_prop:.3e}")
    print(f"\nResulting force ratio (Sun's force / Earth's force): {force_ratio:.2f}")

    print("\nConclusion:")
    if force_ratio > 1:
        print(f"The Sun's gravitational pull on the Moon is about {force_ratio:.2f} times stronger than the Earth's.")
        print("Because the Sun's force is always dominant, the Moon's path is always concave towards the Sun.")
        print("It follows the Earth's orbit while making small 'wobbles', but it never loops backward.")
        print("This corresponds to diagram C.")
    else:
        print("The Earth's pull is stronger, which would lead to loops in the orbit (like diagram D). This is not the case.")

analyze_moon_orbit()
<<<C>>>