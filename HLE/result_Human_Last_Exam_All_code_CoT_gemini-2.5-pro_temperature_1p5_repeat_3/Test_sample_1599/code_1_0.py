import math

def solve_acceleration():
    """
    This function calculates the man's final acceleration based on the problem description.
    It follows a step-by-step process, solving for the parameters of the planned journey first,
    then using that information to model the actual final leg of the journey.
    """

    # --- Initial Given Values ---
    u_initial = 1.5  # m/s
    t0 = 0.0         # s
    t1 = 4.0         # s
    a1 = -0.15       # m/s^2
    t3 = 15.0        # s
    a2 = 0.25        # m/s^2
    gamma_deg = 40.0 # degrees
    t5 = 23.0        # s
    t6 = 40.0        # s
    
    # --- Assumption 1: Resolving the bird's angle ambiguity ---
    # The problem states the bird's ground projection angle is 130 degrees with the man's path (North)
    # but also that the bird moves "northward" and "eastward". This is a contradiction.
    # We assume the intended angle was 180 - 130 = 50 degrees, which allows for North-East motion.
    alpha_deg = 130.0
    beta_deg = 180.0 - alpha_deg # Corrected angle in degrees
    
    # Convert angles to radians for calculations
    beta_rad = math.radians(beta_deg)
    gamma_rad = math.radians(gamma_deg)

    # --- Part 1: Calculate Man's Kinematics for the Planned Journey ---
    
    # Man's state at t1
    y_man_t1 = u_initial * (t1 - t0)
    u_man_t1 = u_initial
    
    # Time t2 is when the man comes to a stop
    t2 = t1 - u_man_t1 / a1
    
    # Man's state at t2 and t3
    y_man_t2 = y_man_t1 + u_man_t1 * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    u_man_t2 = 0.0
    y_man_t3 = y_man_t2 # Man is still from t2 to t3
    u_man_t3 = 0.0
    
    # --- Assumption 2: Man's final planned speed ---
    # We assume the man's speed after accelerating from t3 to t4 returns to his initial speed.
    u_man_t4 = u_initial
    
    # With this, we can find t4 and the man's position at t4
    delta_t_34 = (u_man_t4 - u_man_t3) / a2
    t4 = t3 + delta_t_34
    y_man_t4 = y_man_t3 + u_man_t3 * delta_t_34 + 0.5 * a2 * delta_t_34**2
    
    # --- Part 2: Solve for Bird's Parameters from Planned Rendezvous ---
    # We set up a system of equations for the planned rendezvous where the bird's final
    # x and z coordinates are 0. This allows us to find the bird's speed 'v'.
    # This involves solving a quadratic equation for R = v_g1 / v.
    # The equation is of the form: A*R^2 + B*R + C = 0
    
    C1 = 10 * math.tan(gamma_rad) - (t4 - t3)
    C2 = 4 * math.sin(beta_rad) * math.tan(gamma_rad)
    
    A_quad = 16 + C2**2
    B_quad = 2 * C1 * C2
    C_quad = C1**2 - 16
    
    # Solve quadratic equation for R. We take the positive root.
    R = (-B_quad + math.sqrt(B_quad**2 - 4 * A_quad * C_quad)) / (2 * A_quad)
    
    # Now find the planned final leg duration (Tf) and the bird's speed (v)
    Tf = (10 + 4 * R * math.sin(beta_rad)) / math.cos(gamma_rad)
    v = (y_man_t4 + u_man_t4 * Tf) / (4 * R * math.cos(beta_rad) + (t3-t2))

    v_g1 = R * v # Bird's ground speed in leg 1
    
    # --- Part 3: Determine State of Man and Bird at t5 ---
    
    # Man's state at t5
    # From t4 to t5, man moves at constant speed u_man_t4
    y_man_t5 = y_man_t4 + u_man_t4 * (t5 - t4)
    u_man_t5 = u_man_t4

    # Bird's state at t5
    # First, find bird's position at t4
    x_bird_t4 = 4 * v_g1 * math.sin(beta_rad) + v * (t2 - t1)
    y_bird_t4 = 4 * v_g1 * math.cos(beta_rad) + v * (t3 - t2)
    vz_1 = math.sqrt(v**2 - v_g1**2)
    z_bird_t4 = 4 * vz_1 + v * (t4-t3)
    
    # Then, calculate position at t5, which is (t5-t4) seconds into the planned final leg
    vx_leg5 = -v * math.cos(gamma_rad)
    vz_leg5 = -v * math.sin(gamma_rad)
    
    x_bird_t5 = x_bird_t4 + vx_leg5 * (t5 - t4)
    y_bird_t5 = y_bird_t4 # No planned y-motion in this leg
    z_bird_t5 = z_bird_t4 + vz_leg5 * (t5 - t4)

    # --- Part 4: Model the Actual Final Leg (t5 to t6) ---
    
    delta_t_56 = t6 - t5
    
    # The bird's new velocity components (v_ax, v_ay, v_az)
    # must take it from its position at t5 to the final meeting point (0, y_final, 0)
    v_ax = -x_bird_t5 / delta_t_56
    v_az = -z_bird_t5 / delta_t_56
    # The bird's speed magnitude 'v' is constant
    v_ay_sq = v**2 - v_ax**2 - v_az**2
    v_ay = math.sqrt(v_ay_sq) # Northward component is positive

    # Determine the final meeting point on the y-axis
    y_final = y_bird_t5 + v_ay * delta_t_56

    # --- Part 5: Calculate the Man's Final Acceleration (a3) ---
    # Use the kinematic equation: y_final = y_man_t5 + u_man_t5*t + 0.5*a3*t^2
    # We solve for a3.
    
    displacement = y_final - y_man_t5
    a3 = (displacement - u_man_t5 * delta_t_56) / (0.5 * delta_t_56**2)

    # --- Final Output ---
    print("To find the man's final acceleration, we use the kinematic equation:")
    print(f"y_final = y_initial + v_initial * t + 0.5 * a3 * t^2")
    print("\nBased on the trajectory calculations, the values for this equation are:")
    print(f"Final position (y_final): {y_final:.4f} m")
    print(f"Man's initial position at t5 (y_initial): {y_man_t5:.4f} m")
    print(f"Man's initial velocity at t5 (v_initial): {u_man_t5:.4f} m/s")
    print(f"Time duration (t): {delta_t_56:.4f} s")
    print("\nSolving for a3:")
    print(f"({y_final:.4f}) = ({y_man_t5:.4f}) + ({u_man_t5:.4f}) * ({delta_t_56:.4f}) + 0.5 * a3 * ({delta_t_56:.4f})^2")
    print(f"The required constant acceleration 'a3' is {a3:.4f} m/s^2.")
    
    return a3

# Execute the function and store the final answer
final_acceleration = solve_acceleration()

# The final answer in the required format
print(f"\n<<< {final_acceleration:.4f} >>>")
