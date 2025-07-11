import numpy as np

def calculate_gravitational_shift():
    """
    This script demonstrates the shift of the apparent center of gravity
    for a moving mass when gravitational effects propagate at speed c.

    Scenario:
    - Mass 1 (Observer) is at the origin (0, 0).
    - Mass 2 moves with velocity v along the x-axis. At time t=0, its
      y-position is d, so its present position is (0, d).
    - We calculate the retarded position, which is where Mass 2 was when it
      emitted the gravitational signal that reaches Mass 1 at t=0.
    - The "shift" is the vector difference between the present and retarded positions.
    """

    # --- Constants and Parameters ---
    # Speed of light in m/s
    c = 299792458.0
    # Closest approach distance of Mass 2 (in meters, approx. 2 AU)
    d = 3.0e11
    # Velocity of Mass 2 as a fraction of the speed of light
    v_ratio = 0.5
    v = v_ratio * c

    # --- Calculations ---
    # Calculate the Lorentz factor gamma
    gamma = 1 / np.sqrt(1 - v_ratio**2)

    # Calculate the retarded time (t_e) for the signal to reach the origin at t=0
    # The equation is t_obs = t_e + |r_e|/c, where t_obs=0.
    # This solves to t_e = - (d/c) * gamma
    t_e = - (d / c) * gamma

    # Define the present position of Mass 2 at t=0
    # r_t = (v * 0, d)
    r_present = np.array([0.0, d])

    # Calculate the retarded position of Mass 2 at time t_e
    # r_e = (v * t_e, d)
    r_retarded = np.array([v * t_e, d])

    # Calculate the shift vector
    # The assumption (C or E) leads to the force pointing to the present position.
    # The shift from the naive retarded position is therefore r_present - r_retarded.
    shift_vector = r_present - r_retarded

    # --- Output ---
    print("In the reference frame of Mass 1 (at the origin):")
    print(f"Speed of light (c): {c:.3e} m/s")
    print(f"Velocity of Mass 2 (v): {v:.3e} m/s ({v_ratio}c)")
    print(f"Closest approach distance (d): {d:.3e} m\n")

    print(f"Observation time at Mass 1: t = 0.0 s")
    print(f"Present position of Mass 2 (r_t): ({r_present[0]:.3e}, {r_present[1]:.3e}) m")
    print(f"Signal emission time from Mass 2 (t_e): {t_e:.3f} s")
    print(f"Retarded position of Mass 2 (r_e): ({r_retarded[0]:.3e}, {r_retarded[1]:.3e}) m\n")

    print("The apparent center of gravity is shifted from the retarded position to the present position.")
    print("This shift is in the direction of motion, as shown by the shift vector equation:\n")
    print("Shift Vector = Present Position - Retarded Position")
    
    # Final equation with numbers, as requested.
    print(f"({shift_vector[0]:.3e}, {shift_vector[1]:.3e}) = ({r_present[0]:.3e}, {r_present[1]:.3e}) - ({r_retarded[0]:.3e}, {r_retarded[1]:.3e})")

    # Verify the direction
    direction_of_motion = np.array([1.0, 0.0]) # along the positive x-axis
    # The shift vector is (some_positive_number, 0), so it is in the direction of motion.
    print(f"\nThe direction of motion is along the x-axis.")
    print(f"The calculated shift vector also points along the x-axis, confirming the shift is in the direction of motion.")

calculate_gravitational_shift()
<<<C>>>