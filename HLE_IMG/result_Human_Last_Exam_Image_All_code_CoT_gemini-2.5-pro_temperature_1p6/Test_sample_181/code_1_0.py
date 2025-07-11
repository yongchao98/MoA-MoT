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
    # M_moon is not needed as it cancels out in the ratio
    
    # Average distances
    R_earth_sun = 1.496e11 # Average distance from Earth to Sun (1 AU) in meters
    R_moon_earth = 3.844e8  # Average distance from Moon to Earth in meters

    # The formulas for the forces are:
    # F_sun_on_moon = G * M_sun * M_moon / R_sun_moon^2
    # F_earth_on_moon = G * M_earth * M_moon / R_moon_earth^2
    # The G and M_moon terms cancel when taking the ratio.
    # We approximate the Sun-Moon distance with the Earth-Sun distance.
    
    force_ratio = (M_sun / M_earth) * (R_moon_earth / R_earth_sun)**2
    
    print("Step-by-step reasoning:")
    print("1. We compare the gravitational force from the Sun on the Moon (F_sun) with the force from the Earth on the Moon (F_earth).")
    print(f"2. The ratio of forces F_sun / F_earth is calculated by (M_sun / M_earth) * (R_moon_earth / R_earth_sun)^2.")
    print(f"3. Plugging in the values:")
    print(f"   - Mass of Sun (M_sun) = {M_sun:.3e} kg")
    print(f"   - Mass of Earth (M_earth) = {M_earth:.3e} kg")
    print(f"   - Earth-Sun distance (R_earth_sun) = {R_earth_sun:.3e} m")
    print(f"   - Moon-Earth distance (R_moon_earth) = {R_moon_earth:.3e} m")
    print("")
    print("Calculation:")
    print(f"Force Ratio = ({M_sun:.3e} / {M_earth:.3e}) * ({R_moon_earth:.3e} / {R_earth_sun:.3e})^2")
    print(f"Force Ratio = {M_sun/M_earth:.3f} * ({R_moon_earth/R_earth_sun:.3e})^2")
    print(f"Force Ratio = {M_sun/M_earth:.3f} * {(R_moon_earth/R_earth_sun)**2:.3e}")
    print(f"Final Force Ratio F_sun / F_earth ≈ {force_ratio:.2f}")
    print("\nConclusion:")
    print(f"The Sun's gravitational pull on the Moon is about {force_ratio:.2f} times stronger than the Earth's.")
    print("This means the Moon's path is primarily an orbit around the Sun, which is slightly perturbed by the Earth.")
    print("The trajectory is a 'wavy' line that is always curved (concave) towards the Sun. It never loops or forms sharp points.")
    print("This corresponds to diagram C.")

analyze_moon_orbit()