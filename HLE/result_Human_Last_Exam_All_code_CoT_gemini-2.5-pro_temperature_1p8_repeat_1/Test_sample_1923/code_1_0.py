import numpy as np

def calculate_shift(distance_of_closest_approach, velocity_y, speed_of_gravity):
    """
    Calculates the retarded position and the shift vector for an object
    in transverse motion, as seen by an observer at the origin.

    Args:
        distance_of_closest_approach (float): The distance 'd' of the object at t=0.
                                             Instantaneous position is (d, 0).
        velocity_y (float): The velocity of the object in the y-direction.
                           Velocity vector is (0, v).
        speed_of_gravity (float): The propagation speed of gravity, 'c'.
    """
    if np.abs(velocity_y) >= speed_of_gravity:
        print("Velocity must be less than the speed of gravity.")
        return

    # Instantaneous position at t=0
    p_now = np.array([distance_of_closest_approach, 0.0])

    # Velocity vector
    v = np.array([0.0, velocity_y])

    # Relativistic factors
    beta = velocity_y / speed_of_gravity
    gamma = 1 / np.sqrt(1 - beta**2)

    # The time t_ret when the signal was emitted to be received at t=0
    # For transverse motion, t_ret = -gamma * d / c
    t_ret = -gamma * distance_of_closest_approach / speed_of_gravity

    # Calculate the retarded position: P_ret = P(t=0) + v * t_ret
    # Since P_now is position at t=0, P_ret = P_now + v*t_ret
    p_ret = p_now + v * t_ret

    # The shift vector is the difference between the GR position (P_now)
    # and the naive retarded position (P_ret)
    shift_vector = p_now - p_ret

    print("--- Simulation Setup ---")
    print(f"Observer at: (0, 0)")
    print(f"Object velocity v: {v} m/s (in y-direction)")
    print(f"Speed of gravity c: {speed_of_gravity} m/s")
    print(f"Object's instantaneous position at t=0 (P_now): {p_now} m")
    print("\n--- Calculation Results ---")
    print(f"Signal emitted at retarded time t_ret: {t_ret:.4g} s")
    print(f"Object's retarded position (P_ret): [{p_ret[0]:.4g}, {p_ret[1]:.4g}] m")
    print("\n--- Conclusion ---")
    print("A naive model of finite-speed gravity would point the force towards the Retarded Position (P_ret).")
    print("Assumption (B), describing a consistent relativistic theory like GR, predicts the force points towards the Instantaneous Position (P_now).")
    print(f"The resulting shift in the force's origin is: Shift Vector = P_now - P_ret = [{shift_vector[0]:.4g}, {shift_vector[1]:.4g}]")
    print(f"The object's velocity vector is proportional to {v / np.linalg.norm(v)}.")
    print(f"The shift vector is proportional to {shift_vector / np.linalg.norm(shift_vector)}.")
    # Normalize for comparison
    if np.linalg.norm(shift_vector) > 1e-9:
        is_parallel = np.allclose(v / np.linalg.norm(v), shift_vector / np.linalg.norm(shift_vector))
        if is_parallel:
            print("\nThe shift vector is in the same direction as the velocity vector.")
        else:
            print("\nThe shift vector is NOT in the same direction as the velocity vector.")

# Parameters for the simulation (using realistic values for demonstration)
# Distance to a hypothetical passing object, e.g., 1 AU
d = 1.5e11  # meters
# A relativistic velocity, e.g., 50% the speed of light
v_y = 1.5e8   # m/s
c = 3.0e8     # m/s

calculate_shift(d, v_y, c)