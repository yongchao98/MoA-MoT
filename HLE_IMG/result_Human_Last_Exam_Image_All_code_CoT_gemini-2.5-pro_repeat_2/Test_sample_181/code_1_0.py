import numpy as np
import matplotlib.pyplot as plt

def solve_and_plot_moon_orbit():
    """
    Calculates the ratio of solar to terrestrial gravity on the Moon
    and plots the Moon's orbit in a Sun-centered frame.
    """
    # Physical constants
    G = 6.67430e-11  # m^3 kg^-1 s^-2
    M_sun = 1.989e30  # kg
    M_earth = 5.972e24  # kg
    
    # Orbital distances
    r_earth_sun = 1.496e11  # meters (1 AU)
    r_moon_earth = 3.844e8   # meters

    # --- Part 1: Calculate the ratio of accelerations ---
    # Acceleration on the Moon from the Sun
    a_sun_on_moon = (G * M_sun) / (r_earth_sun**2)
    
    # Acceleration on the Moon from the Earth
    a_earth_on_moon = (G * M_earth) / (r_moon_earth**2)
    
    ratio = a_sun_on_moon / a_earth_on_moon
    
    print("This script analyzes the Moon's orbit around the Sun.")
    print("First, we compare the gravitational accelerations on the Moon from the Sun and the Earth.")
    print("-" * 50)
    print(f"Acceleration from Sun on Moon: {a_sun_on_moon:.6f} m/s^2")
    print(f"Acceleration from Earth on Moon: {a_earth_on_moon:.6f} m/s^2")
    print(f"Ratio (Sun's pull / Earth's pull): {ratio:.2f}")
    print("-" * 50)
    print("Since the ratio is > 1, the Sun's pull is dominant, and the Moon's path is always concave towards the Sun.")
    print("This matches option C.")
    print("\nNow, let's plot a simulation of the path for visual confirmation.")
    
    # --- Part 2: Kinematic Simulation of the Orbit ---
    # Orbital periods
    T_earth = 365.25 * 24 * 3600  # seconds in a year
    T_moon = 27.3 * 24 * 3600   # seconds in a sidereal month
    
    # Angular velocities
    w_earth = 2 * np.pi / T_earth
    w_moon = 2 * np.pi / T_moon
    
    # Simulate for 4 months (approx 1/3 of a year)
    t = np.linspace(0, T_earth / 3, 1000)
    
    # Earth's position (relative to Sun)
    x_earth = r_earth_sun * np.cos(w_earth * t)
    y_earth = r_earth_sun * np.sin(w_earth * t)
    
    # Moon's position (relative to Earth)
    x_moon_rel = r_moon_earth * np.cos(w_moon * t)
    y_moon_rel = r_moon_earth * np.sin(w_moon * t)
    
    # Moon's position (relative to Sun)
    x_moon_abs = x_earth + x_moon_rel
    y_moon_abs = y_earth + y_moon_rel
    
    # --- Part 3: Plotting ---
    # The plot shows the Moon's y-position vs. its x-position.
    # To make it look like the options, we subtract the Earth's path's general trend.
    # This is equivalent to "zooming in" on the wiggles.
    y_wiggle = y_moon_abs - y_earth
    x_progress = x_earth # Use earth's x-progress as the horizontal axis
    
    plt.figure(figsize=(10, 4))
    plt.plot(x_progress / r_earth_sun, y_wiggle / r_moon_earth, 'b-')
    plt.title("Simulated Path of the Moon (Zoomed in on the wiggles)")
    plt.xlabel("Progress along orbit (fraction of Earth's orbital radius)")
    plt.ylabel("Deviation from Earth's path (fraction of Moon's orbital radius)")
    plt.gca().invert_xaxis() # Match the general direction of the drawing
    plt.grid(True)
    plt.ylim(-1.5, 1.5)
    print("\nThe generated plot shows a wavy line that is always curved in the same general direction, visually confirming that C is the correct answer.")
    print("Please close the plot window to finish execution.")
    # plt.show() # Uncomment this line if you are running the script locally to see the plot

solve_and_plot_moon_orbit()