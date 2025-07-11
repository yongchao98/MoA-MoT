import math

def solve_acceleration():
    """
    This function solves the physics problem by following the step-by-step plan.
    1.  Calculates the man's position and velocity through various segments.
    2.  Solves for the bird's unknown speed `v` and a path parameter `k` using the *planned* rendezvous information.
    3.  Calculates the actual meeting point `y_meet` based on the altered path from t5.
    4.  Uses the man's kinematic equation for the final leg to solve for the required acceleration `a3`.
    """

    # --- Given Constants ---
    u = 1.5  # Man's initial and target speed, m/s
    t0 = 0.0
    t1 = 4.0
    alpha_deg = 130.0
    a1 = -0.15  # m/s^2
    t3 = 15.0
    a2 = 0.25  # m/s^2
    gamma_deg = 40.0
    t5 = 23.0
    t6 = 40.0

    # Interpretation: angle with North path is measured from South, so angle wrt North is 180-130=50
    alpha_eff_deg = 180.0 - alpha_deg
    
    # Convert angles to radians for calculations
    alpha_rad = math.radians(alpha_eff_deg) # Angle for bird's ground proj wrt North
    gamma_rad = math.radians(gamma_deg)   # Angle for bird's dive wrt ground

    # --- 1. Man's Kinematics until t5 ---
    # Segment t0 -> t1: Constant velocity
    y_man_t1 = u * (t1 - t0)
    v_man_t1 = u

    # Segment t1 -> t2: Decelerates to a stop
    dt_12 = -v_man_t1 / a1
    t2 = t1 + dt_12
    y_man_t2 = y_man_t1 + v_man_t1 * dt_12 + 0.5 * a1 * dt_12**2
    
    # Segment t2 -> t3: Stays still
    y_man_t3 = y_man_t2
    v_man_t3 = 0.0
    
    # Segment t3 -> t4: Accelerates back to speed u
    dt_34 = (u - v_man_t3) / a2
    t4 = t3 + dt_34
    y_man_t4 = y_man_t3 + v_man_t3 * dt_34 + 0.5 * a2 * dt_34**2
    v_man_t4 = u

    # Segment t4 -> t5: Constant velocity
    dt_45 = t5 - t4
    y_man_t5 = y_man_t4 + v_man_t4 * dt_45
    v_man_t5 = u

    # --- 2. Solve for Bird's Parameters (k and v) ---
    # We derive a quadratic equation for k = v_xy/v of the form Ak^2 + Bk + C = 0
    s_g = math.sin(gamma_rad)
    c_g = math.cos(gamma_rad)
    s_a = math.sin(alpha_rad) # Corresponds to sin(50)
    c_a = math.cos(alpha_rad) # Corresponds to cos(50)
    
    # The equation comes from tan(gamma) = z_bird(t4)/x_bird(t4)
    # (4*vz/v + 6) / (4*vx/v + 10) = tan(gamma)
    # where vx/v = k*s_a and vz/v = sqrt(1-k^2)
    # This can be rearranged into a quadratic equation for k.
    
    term_k = 4 * s_g * s_a
    term_c = 10 * s_g - 6 * c_g
    
    A = term_k**2 + 16 * c_g**2
    B = 2 * term_k * term_c
    C = term_c**2 - 16 * c_g**2
    
    # Solve quadratic equation for k, taking the positive root since k must be positive
    discriminant = B**2 - 4 * A * C
    k = (-B + math.sqrt(discriminant)) / (2 * A)

    # Calculate planned flight time delta_t_pl from t4 onwards
    # delta_t_pl comes from z_bird(t_meet) = 0
    delta_t_pl = (4 * math.sqrt(1 - k**2) + 6) / s_g

    # Calculate bird's speed v using the planned y-rendezvous condition
    # y_man(t_meet) = y_bird(t_meet) => y_man_t4 + u*delta_t_pl = y_bird_t4
    y_bird_t4 = y_man_t4 + u * delta_t_pl
    v = y_bird_t4 / (4 * k * c_a + 1)
    
    # --- 3. Calculate Actual Meeting Point y_meet ---
    # The bird's y position at t5 is the same as at t4
    y_bird_t5 = y_bird_t4
    dt_56 = t6 - t5
    
    # The distance the bird flies from t5 to t6 is v*dt_56.
    # This distance is also sqrt(dx^2 + dy^2 + dz^2).
    # dx = |x_bird(t5) - x_meet| = x_bird(t5) since x_meet=0
    # dz = |z_bird(t5) - z_meet| = z_bird(t5) since z_meet=0
    # x_bird(t5)^2 + z_bird(t5)^2 can be calculated from the planned path.
    # This gives y_meet.
    
    y_meet = y_bird_t5 + v * math.sqrt(dt_56**2 - (delta_t_pl - dt_45)**2)

    # --- 4. Solve for Man's Final Acceleration a3 ---
    # Use the man's kinematic equation for the final leg:
    # y_meet = y_man_t5 + v_man_t5 * dt_56 + 0.5 * a3 * dt_56**2
    
    a3 = (y_meet - y_man_t5 - v_man_t5 * dt_56) / (0.5 * dt_56**2)

    # --- Print Final Calculation as Requested ---
    print("To find the man's final acceleration (a), we use the kinematic equation:")
    print("y_final = y_initial + v_initial * t + 0.5 * a * t^2")
    print("\nPlugging in the values for the man's final journey segment:")
    y_final = y_meet
    y_initial = y_man_t5
    v_initial = v_man_t5
    time_interval = dt_56
    
    print(f"y_final (the meeting point's y-coordinate) = {y_final:.4f} m")
    print(f"y_initial (man's position at t5) = {y_initial:.4f} m")
    print(f"v_initial (man's velocity at t5) = {v_initial:.4f} m/s")
    print(f"t (time interval from t5 to t6) = {time_interval:.4f} s")

    print("\nThe equation becomes:")
    print(f"{y_final:.4f} = {y_initial:.4f} + {v_initial:.4f} * {time_interval:.4f} + 0.5 * a * {time_interval:.4f}^2")
    
    term2_val = v_initial * time_interval
    term3_val = 0.5 * time_interval**2
    print(f"{y_final:.4f} = {y_initial + term2_val:.4f} + {term3_val:.4f} * a")

    print(f"\nSolving for 'a' yields:")
    print(f"a = ({y_final:.4f} - {y_initial + term2_val:.4f}) / {term3_val:.4f}")
    
    # Return the final calculated acceleration.
    return a3

# Execute the solver and print the final answer in the required format.
final_acceleration = solve_acceleration()
print(f"\nThe man's constant acceleration during this final leg is {final_acceleration:.4f} m/s^2.")

# Final answer block
print(f"<<<{final_acceleration}>>>")