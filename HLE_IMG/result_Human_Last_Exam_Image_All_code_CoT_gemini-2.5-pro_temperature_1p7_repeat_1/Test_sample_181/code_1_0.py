import math

def solve_moon_orbit():
    """
    Calculates the ratio of the Sun's gravitational force on the Moon
    to the Earth's gravitational force on the Moon to determine the
    shape of the Moon's orbit in an inertial frame centered on the Sun.
    """

    # --- Physical Constants ---
    # Mass of the Sun in kg
    M_sun = 1.989e30
    # Mass of the Earth in kg
    M_earth = 5.972e24
    # Average distance from Earth to Sun (1 AU) in meters.
    # This is also the approximate distance from the Moon to the Sun.
    r_earth_sun = 1.496e11
    # Average distance from Moon to Earth in meters
    r_moon_earth = 3.844e8

    # --- Calculation ---
    # Gravitational force F = G * m1 * m2 / r^2.
    # We are comparing F_sun_on_moon with F_earth_on_moon.
    # The gravitational constant (G) and the Moon's mass cancel out in the ratio.
    # Ratio = (M_sun / r_earth_sun^2) / (M_earth / r_moon_earth^2)
    # Ratio = (M_sun / M_earth) * (r_moon_earth / r_earth_sun)^2

    mass_ratio = M_sun / M_earth
    distance_ratio_squared = (r_moon_earth / r_earth_sun)**2
    force_ratio = mass_ratio * distance_ratio_squared

    # --- Output ---
    print("Step 1: Define the equation for the ratio of gravitational forces.")
    print("Force Ratio = (Sun's force on Moon) / (Earth's force on Moon)")
    print("Force Ratio = (M_sun / M_earth) * (r_moon_earth / r_earth_sun)^2\n")

    print("Step 2: Substitute the values into the equation.")
    print("M_sun (Mass of Sun) =", f"{M_sun:.3e}", "kg")
    print("M_earth (Mass of Earth) =", f"{M_earth:.3e}", "kg")
    print("r_moon_earth (Distance from Moon to Earth) =", f"{r_moon_earth:.3e}", "m")
    print("r_earth_sun (Distance from Earth to Sun) =", f"{r_earth_sun:.3e}", "m\n")

    print("Step 3: Calculate the final ratio.")
    print("Force Ratio = (", f"{M_sun:.3e}", "/", f"{M_earth:.3e}", ") * (", f"{r_moon_earth:.3e}", "/", f"{r_earth_sun:.3e}", ")^2")
    print(f"Force Ratio = {mass_ratio:.3f} * {distance_ratio_squared:.3e}")
    print(f"Force Ratio ≈ {force_ratio:.2f}\n")

    print("--- Analysis ---")
    print(f"The Sun's pull on the Moon is about {force_ratio:.2f} times stronger than the Earth's pull.")
    print("Because the Sun's gravity is always the dominant force, the Moon's trajectory is always concave towards the Sun.")
    print("This means the path never loops or bends away from the Sun, ruling out options D and E.")
    print("Option C correctly shows a path that is always curved towards the Sun (assumed to be far below the graph), with slight wobbles due to its orbit around the Earth.")

solve_moon_orbit()