import math

def solve_acceleration():
    """
    Calculates the man's final acceleration based on the described scenario.
    """
    # --- Constants ---
    u_initial = 1.5  # m/s
    t0, t1, t3 = 0.0, 4.0, 15.0
    t5, t6 = 23.0, 40.0
    a1 = -0.15  # m/s^2
    a2 = 0.25   # m/s^2
    # The angle alpha=130 from North contradicts "northward".
    # Assuming the angle with the South direction is 130, which means 50 deg from North.
    # Or simply it's a typo for 50 deg. This is required for a physically possible rendezvous.
    alpha_deg = 50.0 # Effective angle from North
    gamma_deg = 40.0 # Angle of bird's dive

    alpha_rad = math.radians(alpha_deg)
    gamma_rad = math.radians(gamma_deg)

    print("Step 1: Calculate the man's position and velocity at key times.")
    # At t1 = 4s
    y_m_t1 = u_initial * t1
    v_m_t1 = u_initial
    print(f"Man's position at t1={t1}s: y = {y_m_t1:.2f} m")

    # Time t2 when man stops
    t2 = t1 + (-v_m_t1 / a1)
    print(f"Man stops decelerating at t2={t2:.2f} s.")
    
    # At t2 = 14s
    y_m_t2 = y_m_t1 + v_m_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    v_m_t2 = 0
    print(f"Man's position at t2={t2:.2f}s: y = {y_m_t2:.2f} m")

    # At t3 = 15s (man is still)
    y_m_t3 = y_m_t2
    v_m_t3 = 0
    print(f"Man's position at t3={t3:.2f}s: y = {y_m_t3:.2f} m")
    
    # Time t4 when man's speed is back to u_initial
    # v_m(t4) = v_m(t3) + a2*(t4-t3) => 1.5 = 0 + 0.25*(t4-15)
    t4 = t3 + u_initial / a2
    print(f"Man reaches {u_initial} m/s again at t4={t4:.2f} s.")
    
    # At t4 = 21s
    y_m_t4 = y_m_t3 + v_m_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
    v_m_t4 = u_initial
    print(f"Man's position at t4={t4:.2f}s: y = {y_m_t4:.2f} m")
    
    print("\nStep 2: Solve for bird's planned flight parameters.")
    # Equation for C_xy (ratio of bird's horizontal speed in seg 1 to total speed)
    # (4*C_xy*sin(alpha) + 10)*tan(gamma) = 4*sqrt(1-C_xy^2) + 6
    # Let's set up the quadratic equation (A^2+C^2)*x^2 + (2AB)*x + (B^2-C^2) = 0 for x=C_xy
    A_quad = 4 * math.sin(alpha_rad) * math.tan(gamma_rad)
    B_quad = 10 * math.tan(gamma_rad) - 6
    C_quad = 4
    
    a_eq = A_quad**2 + C_quad**2
    b_eq = 2 * A_quad * B_quad
    c_eq = B_quad**2 - C_quad**2
    
    # Solve quadratic equation for C_xy
    discriminant = b_eq**2 - 4 * a_eq * c_eq
    C_xy = (-b_eq + math.sqrt(discriminant)) / (2 * a_eq)
    print(f"Solved for bird's speed ratio C_xy = {C_xy:.4f}")

    # Calculate planned meeting time duration delta_t
    # from eq: 4*sqrt(1-C_xy^2) + 6 = delta_t * sin(gamma)
    C_z = math.sqrt(1 - C_xy**2)
    delta_t_meet = (4 * C_z + 6) / math.sin(gamma_rad)
    print(f"Planned meeting time after t4: delta_t_meet = {delta_t_meet:.2f} s")

    # Calculate the planned meeting Y coordinate from man's perspective
    Y_planned = y_m_t4 + v_m_t4 * delta_t_meet
    print(f"Planned meeting location from man's path: Y_planned = {Y_planned:.2f} m")

    # Calculate bird's speed v
    # Y_planned = v * (4 * C_xy * cos(alpha) + 1)
    v_bird = Y_planned / (4 * C_xy * math.cos(alpha_rad) + 1)
    print(f"Bird's constant speed: v = {v_bird:.2f} m/s")

    print("\nStep 3: Calculate Man and Bird positions at t5 = 23s (when plan changes).")
    # Man's position at t5
    y_m_t5 = y_m_t4 + v_m_t4 * (t5 - t4)
    v_m_t5 = v_m_t4
    print(f"Man's state at t5={t5:.0f}s: y = {y_m_t5:.2f} m, v = {v_m_t5:.2f} m/s")

    # Bird's position at t4
    x_b_t4 = v_bird * delta_t_meet * math.cos(gamma_rad)
    y_b_t4 = Y_planned
    z_b_t4 = v_bird * delta_t_meet * math.sin(gamma_rad)
    
    # Bird's position at t5
    time_into_dive = t5 - t4
    x_b_t5 = x_b_t4 - v_bird * math.cos(gamma_rad) * time_into_dive
    y_b_t5 = y_b_t4
    z_b_t5 = z_b_t4 - v_bird * math.sin(gamma_rad) * time_into_dive
    print(f"Bird's position at t5={t5:.0f}s: P_b = ({x_b_t5:.2f}, {y_b_t5:.2f}, {z_b_t5:.2f})")

    print("\nStep 4: Calculate the new meeting point Y_new at t6 = 40s.")
    delta_t_final = t6 - t5
    # Bird must travel from (x_b_t5, y_b_t5, z_b_t5) to (0, Y_new, 0) in delta_t_final
    # Its speed is still v_bird.
    # v_bird^2 = ( (x_b_t5/delta_t_final)^2 + ((Y_new - y_b_t5)/delta_t_final)^2 + (z_b_t5/delta_t_final)^2 )
    # (v_bird * delta_t_final)^2 = x_b_t5^2 + (Y_new - y_b_t5)^2 + z_b_t5^2
    Y_new_minus_y_b_t5_sq = (v_bird * delta_t_final)**2 - x_b_t5**2 - z_b_t5**2
    # Since the bird shifts "northward", Y_new > y_b_t5, so we take the positive root.
    Y_new = y_b_t5 + math.sqrt(Y_new_minus_y_b_t5_sq)
    print(f"New meeting location: Y_new = {Y_new:.2f} m")

    print("\nStep 5: Calculate the man's final acceleration, a3.")
    # Man travels from y_m_t5 to Y_new in delta_t_final with initial velocity v_m_t5
    # delta_y = v0*t + 0.5*a*t^2
    delta_y = Y_new - y_m_t5
    # a3 = (delta_y - v0*t) * 2 / t^2
    a3 = (delta_y - v_m_t5 * delta_t_final) * 2 / (delta_t_final**2)
    print(f"Equation for acceleration: {delta_y:.2f} m = ({v_m_t5:.2f} m/s * {delta_t_final:.2f} s) + 0.5 * a3 * ({delta_t_final:.2f} s)^2")
    print(f"Solving for a3: ({delta_y:.2f} - {v_m_t5 * delta_t_final:.2f}) * 2 / {delta_t_final**2:.2f} = a3")
    print(f"Man's constant acceleration a3 = {a3:.4f} m/s^2")
    
    return a3

final_acceleration = solve_acceleration()
print(f"\nThe man's constant acceleration during the final leg of his journey is {final_acceleration:.4f} m/s^2.")
print(f"<<<{final_acceleration:.4f}>>>")
