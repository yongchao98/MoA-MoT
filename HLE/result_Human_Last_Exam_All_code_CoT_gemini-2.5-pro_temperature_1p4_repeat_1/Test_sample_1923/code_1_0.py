import math

def calculate_gravity_shift():
    """
    Demonstrates how assumption C leads to a forward shift in the center of gravity.
    """
    # --- Simulation Parameters ---
    # Speed of light
    c = 299792458  # m/s
    # Velocity of mass 2 (as a fraction of c)
    v_fraction = 0.5
    v = v_fraction * c  # m/s
    # Perpendicular distance from mass 1 to the path of mass 2's center
    y0 = 1.0e9  # meters (1 million km)
    # Half-length of mass 2 (from its center to its front/back)
    L = 1.0e6  # meters (1,000 km)

    print("--- Simulation Parameters ---")
    print(f"Speed of light (c): {c} m/s")
    print(f"Velocity of mass 2 (v): {v:.2e} m/s ({v_fraction*100}% of c)")
    print(f"Closest approach distance (y0): {y0:.2e} m")
    print(f"Object half-length (L): {L:.2e} m\n")

    # We analyze the moment of closest approach (transverse motion).
    # Observer (mass 1) is at the origin (0, 0).
    # Center of mass 2 is at (0, y0). Velocity is (v, 0).
    
    # Position of the front of mass 2
    pos_front = (L, y0)
    # Position of the back of mass 2
    pos_back = (-L, y0)

    # --- Calculations for the Front ---
    # Vector from the front to the observer
    r_vec_front = (-pos_front[0], -pos_front[1])
    # Current distance from the front to the observer
    d_curr_front = math.sqrt(r_vec_front[0]**2 + r_vec_front[1]**2)
    # Radial velocity of the front (v . r / |r|)
    # v_vec is (v, 0)
    v_rad_front = (v * r_vec_front[0]) / d_curr_front
    # Propagation time from the front
    dt_front = d_curr_front / (c - v_rad_front)
    # Field strength from the front (Strength ~ 1/Δt)
    strength_front = 1 / dt_front

    # --- Calculations for the Back ---
    # Vector from the back to the observer
    r_vec_back = (-pos_back[0], -pos_back[1])
    # Current distance from the back to the observer
    d_curr_back = math.sqrt(r_vec_back[0]**2 + r_vec_back[1]**2)
    # Radial velocity of the back
    v_rad_back = (v * r_vec_back[0]) / d_curr_back
    # Propagation time from the back
    dt_back = d_curr_back / (c - v_rad_back)
    # Field strength from the back (Strength ~ 1/Δt)
    strength_back = 1 / dt_back
    
    print("--- Analysis Results ---")
    print(f"Front of Object:")
    print(f"  - Radial velocity: {v_rad_front:.3e} m/s (negative = approaching)")
    print(f"  - Propagation time (Δt): {dt_front:.4f} s")
    print(f"  - Resulting field strength (1/Δt): {strength_front:.4e}\n")

    print(f"Back of Object:")
    print(f"  - Radial velocity: {v_rad_back:.3e} m/s (positive = receding)")
    print(f"  - Propagation time (Δt): {dt_back:.4f} s")
    print(f"  - Resulting field strength (1/Δt): {strength_back:.4e}\n")

    print("--- Conclusion ---")
    print("The field strength from the front is greater than from the back.")
    print("Final Equation:")
    print(f"{strength_front:.4e} > {strength_back:.4e}")
    is_stronger = "True" if strength_front > strength_back else "False"
    print(f"Is Front Strength > Back Strength? {is_stronger}")

calculate_gravity_shift()