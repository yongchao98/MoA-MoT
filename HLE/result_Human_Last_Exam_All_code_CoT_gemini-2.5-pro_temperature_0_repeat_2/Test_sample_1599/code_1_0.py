import math

def solve_acceleration():
    """
    Calculates the man's final acceleration, a_3.
    """
    # Given constants
    u_initial = 1.5  # m/s
    t0 = 0.0
    t1 = 4.0
    alpha_deg = 130.0
    a1 = -0.15  # m/s^2
    t3 = 15.0
    a2 = 0.25  # m/s^2
    gamma_deg = 40.0
    t5 = 23.0
    t6 = 40.0

    # --- Step 1: Analyze the Man's Journey ---
    # Position and velocity at t1
    y_m_t1 = u_initial * t1
    u_m_t1 = u_initial

    # Time t2 when the man stops
    # u_m_t2 = u_m_t1 + a1 * (t2 - t1) => 0 = 1.5 - 0.15 * (t2 - 4)
    t2 = t1 + (-u_m_t1 / a1)

    # Position at t2
    y_m_t2 = y_m_t1 + u_m_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    u_m_t2 = 0.0

    # Position at t3 (man is still)
    y_m_t3 = y_m_t2
    u_m_t3 = 0.0

    # Time t4 when man reaches u_initial again
    # u_m_t4 = u_m_t3 + a2 * (t4 - t3) => 1.5 = 0 + 0.25 * (t4 - 15)
    t4 = t3 + (u_initial - u_m_t3) / a2

    # Position at t4
    y_m_t4 = y_m_t3 + u_m_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
    u_m_t4 = u_initial

    # Position at t5
    y_m_t5 = y_m_t4 + u_m_t4 * (t5 - t4)
    u_m_t5 = u_initial

    # --- Step 2: Find the Bird's Speed v ---
    # As explained in the plan, we interpret alpha as the angle with the South direction
    # to resolve the contradiction in the problem statement.
    alpha = math.radians(alpha_deg)
    gamma = math.radians(gamma_deg)

    # Let k = v_xy / v, where v_xy is the bird's ground speed in the first segment.
    # From the planned rendezvous, we derive a quadratic equation for k^2.
    # (3.064176 k + 2.8495)^2 = (4.767)^2 * (1-k^2)
    # This simplifies to: 32.113 k^2 + 17.463 k - 14.6044 = 0
    a_k = 32.113
    b_k = 17.463
    c_k = -14.6044
    # Solve the quadratic equation for k (must be positive)
    k = (-b_k + math.sqrt(b_k**2 - 4 * a_k * c_k)) / (2 * a_k)

    # Now use the second rendezvous condition to find v
    # 18 + 2.3335 * (4*sqrt(1-k^2) + 6) = v * (2.57115 * k + 1)
    # Note: these constants are derived from the problem's geometric constraints.
    lhs_v = y_m_t4 + (u_initial / math.sin(gamma)) * (4 * math.sqrt(1 - k**2) + 6)
    rhs_v_factor = (4 * (-math.cos(alpha)) * k + 1)
    v = lhs_v / rhs_v_factor

    # --- Step 3 & 4: Determine Bird's Final Position y_b(t6) ---
    # Bird's velocity components in the first segment (t0 to t1)
    v_xy = k * v
    v_bx1 = v_xy * math.sin(alpha)
    v_by1 = v_xy * (-math.cos(alpha)) # Northward
    v_bz1 = v * math.sqrt(1 - k**2)

    # Bird's position at t4
    x_b_t4 = 4 * v_bx1 + v * (t2 - t1)
    y_b_t4 = 4 * v_by1 + v * (t3 - t2)
    z_b_t4 = 4 * v_bz1 + v * (t4 - t3)

    # Bird's position at t5
    x_b_t5 = x_b_t4 - v * math.cos(gamma) * (t5 - t4)
    y_b_t5 = y_b_t4
    z_b_t5 = z_b_t4 - v * math.sin(gamma) * (t5 - t4)

    # Bird's final y-position at t6
    # The distance traveled in the y-direction is sqrt(total_dist^2 - x_dist^2 - z_dist^2)
    # total_dist = v * (t6-t5), x_dist = x_b_t5, z_dist = z_b_t5
    y_dist_56 = math.sqrt((v * (t6 - t5))**2 - x_b_t5**2 - z_b_t5**2)
    y_b_t6 = y_b_t5 + y_dist_56

    # --- Step 5: Solve for a3 ---
    # Man's final position y_m(t6) = y_b(t6)
    # y_m(t6) = y_m_t5 + u_m_t5 * (t6 - t5) + 0.5 * a3 * (t6 - t5)**2
    delta_t_56 = t6 - t5
    a3 = (y_b_t6 - y_m_t5 - u_m_t5 * delta_t_56) / (0.5 * delta_t_56**2)

    # --- Output the final equation and result ---
    print("The final acceleration 'a3' is found by solving the kinematic equation for the man's final leg of the journey:")
    print(f"y_m(t6) = y_m(t5) + u(t5)*(t6-t5) + 0.5*a3*(t6-t5)^2\n")
    print("Where:")
    print(f"The final rendezvous position, y_m(t6), is determined by the bird's final position, y_b(t6) = {y_b_t6:.3f} m.")
    print(f"The man's position at t5, y_m(t5), is {y_m_t5:.3f} m.")
    print(f"The man's velocity at t5, u(t5), is {u_m_t5:.3f} m/s.")
    print(f"The final time interval is (t6-t5) = ({t6} - {t5}) = {delta_t_56:.3f} s.\n")
    
    print("Plugging in the values:")
    print(f"{y_b_t6:.3f} = {y_m_t5:.3f} + {u_m_t5:.3f} * {delta_t_56:.3f} + 0.5 * a3 * ({delta_t_56:.3f})^2")
    print(f"{y_b_t6:.3f} = {y_m_t5 + u_m_t5 * delta_t_56:.3f} + {0.5 * delta_t_56**2:.3f} * a3")
    print(f"{y_b_t6 - (y_m_t5 + u_m_t5 * delta_t_56):.3f} = {0.5 * delta_t_56**2:.3f} * a3")
    print(f"a3 = {y_b_t6 - (y_m_t5 + u_m_t5 * delta_t_56):.3f} / {0.5 * delta_t_56**2:.3f}\n")
    
    print(f"The man's constant acceleration during the final leg is a3 = {a3:.3f} m/s^2.")
    
    return a3

# Run the solver
final_acceleration = solve_acceleration()
print(f"\n<<< {final_acceleration:.3f} >>>")