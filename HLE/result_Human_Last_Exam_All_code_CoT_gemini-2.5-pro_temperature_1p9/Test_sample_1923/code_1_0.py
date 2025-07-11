import math

def solve_gravitational_aberration():
    """
    Calculates the apparent shift of a moving mass due to the finite speed of gravity
    in a naive model where the force points to the retarded position.
    """
    # --- Parameters ---
    c = 299792458.0  # Speed of light in m/s
    v = 0.5 * c      # Velocity of mass 2 (50% speed of light)
    h = 1.0e11       # Perpendicular distance between mass 1's and mass 2's paths (in meters)
    t = 1000.0       # Time of observation (in seconds)

    print(f"--- System Parameters ---")
    print(f"Observer: Mass 1 at (0, 0)")
    print(f"Source: Mass 2 moving at v = {v:.2e} m/s (0.5c)")
    print(f"Path of Mass 2: y = {h:.2e} m")
    print(f"Time of observation: t = {t} s\n")

    # --- Calculations ---
    
    # 1. Instantaneous position of mass 2 at time t
    x_inst = v * t
    y_inst = h
    print(f"--- Positions at t = {t} s ---")
    print(f"Instantaneous Position of Mass 2: ({x_inst:.2e}, {y_inst:.2e}) m")

    # 2. Find the retarded time t' for the signal to reach mass 1 at time t
    # We solve the quadratic equation for t': (c²-v²)t'² - (2c²t)t' + (c²t²-h²) = 0
    a = c**2 - v**2
    b = -2 * c**2 * t
    c_const = c**2 * t**2 - h**2

    # Using the quadratic formula, and taking the solution t' < t
    discriminant = b**2 - 4 * a * c_const
    if discriminant < 0:
        print("Error: No real solution for retarded time. The event is not yet visible.")
        return

    t_retarded = (-b - math.sqrt(discriminant)) / (2 * a)

    # 3. Retarded position of mass 2 at time t'
    x_retarded = v * t_retarded
    y_retarded = h
    print(f"Retarded Position of Mass 2:     ({x_retarded:.2e}, {y_retarded:.2e}) m (as perceived at t={t}s)")

    # 4. Calculate the shift vector and its direction
    shift_x = x_retarded - x_inst
    shift_y = y_retarded - y_inst
    
    print("\n--- Aberration Result ---")
    # Propagation time = t - t_retarded
    propagation_time = t - t_retarded
    print(f"The signal took {propagation_time:.2f} seconds to travel from Mass 2 to Mass 1.")

    print(f"The apparent position of Mass 2 is shifted from its actual position.")
    print(f"Shift Vector: ({shift_x:.2e}, {shift_y:.2e}) m")

    # Check direction of the shift relative to motion
    direction_of_motion = "positive X" if v > 0 else "negative X"
    shift_direction = "positive X" if shift_x > 0 else "negative X"

    print(f"\nDirection of Motion: {direction_of_motion}")
    print(f"Direction of Shift:  {shift_direction}")

    if (v > 0 and shift_x < 0) or (v < 0 and shift_x > 0):
        print("\nConclusion: The center of gravity appears shifted in the direction OPPOSITE to its motion.")
    else:
        print("\nConclusion: The center of gravity appears shifted IN the direction of its motion.")

solve_gravitational_aberration()

<<<C>>>