import math

def solve_moon_orbit_shape():
    """
    Calculates and explains the physical parameters that determine the shape
    of the Moon's orbit in a heliocentric frame.
    """
    # Physical constants
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_sun = 1.989e30  # Mass of the Sun in kg
    M_earth = 5.972e24 # Mass of the Earth in kg
    M_moon = 7.342e22  # Mass of the Moon in kg

    # Average orbital radii
    R_earth_sun = 1.496e11 # Earth-Sun distance (1 AU) in meters
    R_moon_earth = 3.844e8  # Moon-Earth distance in meters

    # --- Step 1: Compare the gravitational forces on the Moon ---
    # The formula for gravitational force is F = G * M1 * M2 / r^2

    # Force exerted by the Sun on the Moon
    F_sun_on_moon = (G * M_sun * M_moon) / R_earth_sun**2

    # Force exerted by the Earth on the Moon
    F_earth_on_moon = (G * M_earth * M_moon) / R_moon_earth**2

    # Calculate the ratio
    force_ratio = F_sun_on_moon / F_earth_on_moon

    print("--- Analysis of the Moon's Orbit ---")
    print("\nStep 1: Comparing Gravitational Forces on the Moon")
    print(f"The Sun's gravitational force on the Moon is calculated as: F_sun_moon = G * M_sun * M_moon / R_sun_moon^2")
    print(f"F_sun_moon = {G:.3e} * {M_sun:.3e} * {M_moon:.3e} / ({R_earth_sun:.3e})^2 = {F_sun_on_moon:.3e} N")
    print(f"The Earth's gravitational force on the Moon is calculated as: F_earth_moon = G * M_earth * M_moon / R_earth_moon^2")
    print(f"F_earth_moon = {G:.3e} * {M_earth:.3e} * {M_moon:.3e} / ({R_moon_earth:.3e})^2 = {F_earth_on_moon:.3e} N")
    print(f"\nThe ratio of these forces (Sun's force / Earth's force) is {F_sun_on_moon:.3e} / {F_earth_on_moon:.3e} = {force_ratio:.2f}")
    print("\nConclusion from forces:")
    print(f"The Sun's gravitational pull on the Moon is about {force_ratio:.2f} times STRONGER than the Earth's pull.")
    print("This means the Moon's path is always accelerating towards the Sun. Therefore, its trajectory is always concave towards the Sun.")
    print("This rules out options D and E, which show the path bending away from the Sun (creating loops or cusps).")

    # --- Step 2: Compare the orbital velocities ---
    # The formula for orbital velocity is v = 2 * pi * r / T
    T_earth_orbit_s = 365.25 * 24 * 3600 # Earth's orbital period in seconds
    T_moon_orbit_s = 27.3 * 24 * 3600    # Moon's orbital period (sidereal) in seconds

    # Earth's orbital velocity around the Sun
    V_earth_around_sun = (2 * math.pi * R_earth_sun) / T_earth_orbit_s

    # Moon's orbital velocity around the Earth
    V_moon_around_earth = (2 * math.pi * R_moon_earth) / T_moon_orbit_s

    print("\nStep 2: Comparing Orbital Velocities")
    print(f"The Earth's orbital velocity around the Sun is v_earth = 2 * pi * {R_earth_sun:.3e} / {T_earth_orbit_s:.3e} = {V_earth_around_sun:,.0f} m/s.")
    print(f"The Moon's orbital velocity around the Earth is v_moon = 2 * pi * {R_moon_earth:.3e} / {T_moon_orbit_s:.3e} = {V_moon_around_earth:,.0f} m/s.")
    print("\nConclusion from velocities:")
    print(f"The Moon is carried along by the Earth at a high speed (~{V_earth_around_sun:,.0f} m/s), while its own orbital motion around the Earth is much slower (~{V_moon_around_earth:,.0f} m/s).")
    print("Because the Moon's speed around the Earth is much less than the Earth's speed around the Sun, the Moon's overall velocity never reverses. It just speeds up and slows down slightly along its path.")
    print("This confirms that the path cannot have loops, as shown in D.")

    # --- Step 3: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    print("1. The Moon orbits the Earth, so its path is not a simple smooth curve (ruling out A and B). It must have a 'wobble'.")
    print("2. The Moon's path is always concave towards the Sun because the Sun's gravity is the dominant force (ruling out D and E).")
    print("\nTherefore, the correct representation is a large-scale orbit that is always concave towards the Sun, with a smaller, periodic wobble superimposed on it.")
    print("This perfectly matches option C.")

if __name__ == '__main__':
    solve_moon_orbit_shape()