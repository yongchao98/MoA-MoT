import math

def solve_acceleration():
    """
    Calculates the man's constant acceleration during the final leg of his journey.
    """
    # Given constants
    u = 1.5  # m/s, man's initial speed
    t0 = 0.0
    t1 = 4.0
    t3 = 15.0
    t5 = 23.0
    t6 = 40.0
    a1 = -0.15  # m/s^2
    a2 = 0.25  # m/s^2
    alpha_deg = 130.0 # degrees
    gamma_deg = 40.0 # degrees

    # Step 1: Analyze Man's motion to find t2 and t4
    # From t1 to t2, man decelerates from u to 0
    # v = u_at_t1 + a1 * (t2 - t1) => 0 = 1.5 + (-0.15) * (t2 - 4)
    t2 = t1 - u / a1
    
    # From t3 to t4, man accelerates from 0 to u
    # v = u_at_t3 + a2 * (t4 - t3) => 1.5 = 0 + 0.25 * (t4 - 15)
    t4 = t3 + u / a2

    # Step 2 & 3: Model bird's motion and solve for its speed 'v'
    # Based on the planned meeting, we derive equations for bird's speed v and flight parameters.
    # A contradiction arises if alpha is measured from North. We assume alpha is measured from East.
    # Let k = v_xy1 / v, where v_xy1 is the bird's horizontal speed in the first segment.
    
    # Trigonometric values
    alpha_rad = math.radians(alpha_deg)
    gamma_rad = math.radians(gamma_deg)
    sin_alpha = math.sin(alpha_rad)
    cos_alpha = math.cos(alpha_rad)
    cos_gamma = math.cos(gamma_rad)
    tan_gamma = math.tan(gamma_rad)
    sin_gamma = math.sin(gamma_rad)

    # From planned meeting conditions, we get a quadratic equation for k:
    # A*k^2 + B*k + C = 0
    # derived from ( (4*k*cos(alpha) + 10)*tan(gamma) - 6 )^2 = 16*(1-k^2)
    
    # Coefficients for the quadratic equation in k
    term1 = 4 * cos_alpha * tan_gamma
    term2 = 10 * tan_gamma - 6
    
    A_k = term1**2 + 16
    B_k = 2 * term1 * term2
    C_k = term2**2 - 16
    
    # Solve for k using the quadratic formula
    discriminant_k = B_k**2 - 4 * A_k * C_k
    k = (-B_k + math.sqrt(discriminant_k)) / (2 * A_k)

    # Now solve for v using the other condition from the planned meeting
    # v * (4*k*sin(alpha) + 1) = 18 + (1.5/cos(gamma)) * (4*k*cos(alpha) + 10)
    v_num = 18 + (u / cos_gamma) * (4 * k * cos_alpha + 10)
    v_den = 4 * k * sin_alpha + 1
    v = v_num / v_den

    # Step 4: Determine positions at t5
    # Man's position at t4
    y_man_1 = u * t1
    y_man_2 = y_man_1 + u * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    y_man_3 = y_man_2
    y_man_4 = y_man_3 + 0.5 * a2 * (t4 - t3)**2
    
    # Man's position and velocity at t5
    y_man_5 = y_man_4 + u * (t5 - t4)
    u_man_5 = u

    # Bird's position at t4
    v_xy1 = k * v
    v_z1 = v * math.sqrt(1 - k**2)
    
    x_bird_4 = 4 * v_xy1 * cos_alpha + 10 * v
    y_bird_4 = 4 * v_xy1 * sin_alpha + v
    z_bird_4 = 4 * v_z1 + 6 * v

    # Bird's position at t5 (2 seconds into the planned dive)
    dt_4_5 = t5 - t4
    x_bird_5 = x_bird_4 - v * cos_gamma * dt_4_5
    y_bird_5 = y_bird_4
    z_bird_5 = z_bird_4 - v * sin_gamma * dt_4_5

    # Step 5 & 6: Analyze final leg and find meeting point
    dt_5_6 = t6 - t5
    # Bird's total distance squared from t5 to t6
    dist_sq_bird = (v * dt_5_6)**2
    
    # (y_man_6 - y_bird_5)^2 = dist_sq_bird - x_bird_5^2 - z_bird_5^2
    y_diff_sq = dist_sq_bird - x_bird_5**2 - z_bird_5**2
    y_diff = math.sqrt(y_diff_sq)
    
    # The meeting point must be north of the bird's y-position
    y_man_6 = y_bird_5 + y_diff
    
    # Step 7: Calculate the man's final acceleration
    # y_man_6 = y_man_5 + u_man_5 * dt_5_6 + 0.5 * a3 * dt_5_6**2
    
    numerator = y_man_6 - y_man_5 - u_man_5 * dt_5_6
    denominator = 0.5 * dt_5_6**2
    a3 = numerator / denominator

    print("Calculation of the man's final acceleration, a3.")
    print("The final acceleration is calculated using the formula:")
    print("a3 = (y_man_6 - y_man_5 - u_man_5 * dt) / (0.5 * dt^2)\n")
    print("Where:")
    print(f"y_man_6 (final meeting y-position) = {y_man_6:.4f} m")
    print(f"y_man_5 (man's y-position at t5) = {y_man_5:.4f} m")
    print(f"u_man_5 (man's speed at t5) = {u_man_5:.4f} m/s")
    print(f"dt (time interval from t5 to t6) = {dt_5_6:.4f} s\n")
    print("Substituting the values into the equation:")
    print(f"a3 = ({y_man_6:.4f} - {y_man_5:.4f} - {u_man_5:.4f} * {dt_5_6:.4f}) / (0.5 * {dt_5_6:.4f}^2)")
    print(f"a3 = ({numerator:.4f}) / ({denominator:.4f})")
    print(f"a3 = {a3:.4f} m/s^2")
    
    return a3

# Execute the function and store the result
final_acceleration = solve_acceleration()

# Final answer format
print(f"\n<<<{final_acceleration:.4f}>>>")