import numpy as np
import matplotlib.pyplot as plt

def solve_and_plot():
    """
    Calculates key physical parameters and plots the Moon's orbit around the Sun.
    """
    # --- Constants (in SI units: meters, kg, seconds) ---
    G = 6.67430e-11  # Gravitational constant
    M_sun = 1.989e30  # Mass of the Sun
    M_earth = 5.972e24 # Mass of the Earth
    M_moon = 7.342e22  # Mass of the Moon

    # Distances (average orbital radii)
    R_earth_sun = 1.496e11 # Earth-Sun distance (1 AU)
    R_moon_earth = 3.844e8  # Moon-Earth distance

    # Orbital Periods
    T_earth = 365.25 * 24 * 3600 # Earth's orbital period (1 year)
    T_moon = 27.3 * 24 * 3600    # Moon's sidereal orbital period

    # --- Step 1: Force Calculation ---
    F_sun_on_moon = G * M_sun * M_moon / R_earth_sun**2
    F_earth_on_moon = G * M_earth * M_moon / R_moon_earth**2
    force_ratio = F_sun_on_moon / F_earth_on_moon

    print("--- Physics Analysis ---")
    print(f"Gravitational force of Sun on Moon: {F_sun_on_moon:.2e} N")
    print(f"Gravitational force of Earth on Moon: {F_earth_on_moon:.2e} N")
    print(f"Ratio (Sun's force / Earth's force): {force_ratio:.2f}")
    print("The Sun's gravitational pull on the Moon is more than twice the Earth's pull.\n")

    # --- Step 2: Velocity Calculation ---
    V_earth = 2 * np.pi * R_earth_sun / T_earth
    V_moon_around_earth = 2 * np.pi * R_moon_earth / T_moon

    print(f"Earth's orbital velocity around Sun: {V_earth/1000:.2f} km/s")
    print(f"Moon's orbital velocity around Earth: {V_moon_around_earth/1000:.2f} km/s")
    print(f"Ratio (Earth's speed / Moon's speed): {V_earth/V_moon_around_earth:.2f}")
    print("The Earth's orbital speed is about 30 times faster than the Moon's.\n")
    print("Conclusion: The Moon's path is always concave towards the Sun and has no loops or cusps.")
    print("The resulting path resembles a wave on a large arc, matching option C.")


    # --- Step 3: Simulation and Plotting ---
    # Time array for 4 months (approx 1/3 of a year)
    t = np.linspace(0, T_earth / 3, 1000)

    # Angular velocities
    omega_earth = 2 * np.pi / T_earth
    omega_moon = 2 * np.pi / T_moon

    # Parametric equations for positions (Sun at origin)
    # Earth's position
    x_earth = R_earth_sun * np.cos(omega_earth * t)
    y_earth = R_earth_sun * np.sin(omega_earth * t)

    # Moon's position relative to Earth
    x_moon_rel = R_moon_earth * np.cos(omega_moon * t)
    y_moon_rel = R_moon_earth * np.sin(omega_moon * t)

    # Moon's absolute position
    x_moon = x_earth + x_moon_rel
    y_moon = y_earth + y_moon_rel

    # Plotting
    plt.figure(figsize=(10, 8))
    plt.plot(x_earth, y_earth, label="Earth's Orbit", linestyle='--', color='blue')
    plt.plot(x_moon, y_moon, label="Moon's Orbit around Sun", color='gray', linewidth=2)
    
    # We need to zoom in significantly to see the wobble
    # Let's adjust the aspect ratio to make the small wobble visible
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    
    # We will plot a much smaller segment and adjust aspect ratio to see the shape
    # Since the arc is so large, plotting a segment and "straightening" it reveals the wave
    # Let's plot y vs x but with a highly distorted aspect ratio to see the wave
    plt.figure(figsize=(12, 4)) # Wide and short figure
    plt.plot(x_moon, y_moon, label="Moon's Path (distorted aspect ratio)", color='black')
    plt.title("Visualisation of the Moon's Path (Not to Scale)")
    plt.xlabel("Position X (m)")
    plt.ylabel("Position Y (m)")
    plt.tick_params(axis='y', labelleft=False, left=False) # Hide y-axis for clarity
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    
    print("\nGenerating plot to visualize the path shape (Note: Aspect ratio is distorted to make the wobble visible)...")
    plt.show()

solve_and_plot()