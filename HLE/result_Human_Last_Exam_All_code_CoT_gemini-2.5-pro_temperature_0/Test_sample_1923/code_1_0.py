import numpy as np

def calculate_positions(v, b, c=1.0):
    """
    Calculates the instantaneous and retarded positions of a moving object.

    The scenario:
    - An observer is at the origin (0, 0).
    - A source object moves with velocity v along the line y = b.
    - We calculate the positions at the moment the object's true position is (0, b).

    Args:
        v (float): The velocity of the source object (as a fraction of c).
        b (float): The distance of closest approach to the observer.
        c (float): The speed of light (and gravity).
    """
    if v >= c:
        print("Error: Velocity must be less than the speed of light.")
        return

    # The instantaneous position of the source at time t=0
    r_inst = np.array([0, b])

    # Calculate the Lorentz factor gamma
    gamma = 1 / np.sqrt(1 - (v/c)**2)

    # The signal that reaches the observer at t=0 was emitted at an earlier, "retarded" time.
    # For this geometry, the retarded time t_ret = -gamma * b / c
    t_ret = -gamma * b / c

    # The retarded position is the position of the source at the retarded time.
    # x_ret = v * t_ret
    r_ret = np.array([v * t_ret, b])

    # The "shift" described in the problem can be interpreted as the vector
    # from the retarded position to the instantaneous position.
    shift_vector = r_inst - r_ret

    print(f"Setup:")
    print(f"  - Velocity of mass 2 (v): {v:.2f}c")
    print(f"  - Distance of closest approach (b): {b:.2f} units")
    print("-" * 30)
    print(f"Results (at observation time t=0):")
    print(f"  - Instantaneous Position of mass 2: ({r_inst[0]:.3f}, {r_inst[1]:.3f})")
    print(f"  - Retarded Position of mass 2:    ({r_ret[0]:.3f}, {r_ret[1]:.3f})")
    print("-" * 30)
    print(f"Analysis of the Shift:")
    print(f"  - The shift vector (Instantaneous - Retarded) is ({shift_vector[0]:.3f}, {shift_vector[1]:.3f})")
    print(f"  - The direction of motion is along the x-axis (e.g., ({v:.2f}, 0.00)).")
    print("\nConclusion:")
    print("The apparent center of gravity (the instantaneous position) is shifted from the naive retarded position.")
    print("This shift vector points purely in the direction of the object's motion.")

# --- Parameters for the demonstration ---
# Velocity of mass 2 as a fraction of the speed of light c
velocity = 0.6 * 1.0 # 60% of the speed of light

# Distance of closest approach
distance_b = 10.0

calculate_positions(v=velocity, b=distance_b)