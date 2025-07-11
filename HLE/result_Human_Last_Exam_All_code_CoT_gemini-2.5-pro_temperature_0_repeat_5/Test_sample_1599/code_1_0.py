import math

def solve_acceleration():
    """
    This function calculates the man's final acceleration based on the described scenario.
    It follows a step-by-step process to determine the bird's speed and the positions
    of both the man and the bird at key moments to find the final answer.
    """
    # Step 1: Define constants from the problem
    u_man_initial = 1.5  # m/s
    t1 = 4.0  # s
    a1 = -0.15  # m/s^2
    t3 = 15.0  # s
    a2 = 0.25  # m/s^2
    gamma_rad = math.radians(40.0)
    t5 = 23.0  # s
    t6 = 40.0  # s
    
    # Angle alpha is interpreted as 130 degrees from the East axis.
    alpha_rad = math.radians(130.0)

    # Step 2: Analyze the man's motion to find t2
    # Man's velocity becomes 0 at t2. u(t2) = u(t1) + a1*(t2-t1) = 0
    # u(t1) is the initial constant speed.
    # 0 = 1.5 + (-0.15) * (t2 - 4) => 1.5 = 0.15 * (t2 - 4) => 10 = t2 - 4
    t2 = 14.0  # s

    # Step 3: Use the planned rendezvous to find the bird's speed 'v'
    # We assume the bird's initial motion is purely horizontal (z-component is zero)
    # to resolve ambiguity in the problem statement.
    
    # Bird's position at t4 in terms of v and t4
    # x_b(t4) = x_b(t1) + v*(t2-t1) = 4*v*cos(alpha) + v*10
    x_b_t4_factor = 4 * math.cos(alpha_rad) + 10
    # y_b(t4) = y_b(t1) + v*(t3-t2) = 4*v*sin(alpha) + v*1
    y_b_t4_factor = 4 * math.sin(alpha_rad) + 1
    
    # From z_b(t4) = x_b(t4) * tan(gamma), we find t4
    # v*(t4-t3) = (v*x_b_t4_factor) * tan(gamma)
    t4 = t3 + x_b_t4_factor * math.tan(gamma_rad)
    
    # Now find the man's position at the planned rendezvous time tr
    # tr - t4 = x_b(t4) / (v*cos(gamma))
    delta_t_r = x_b_t4_factor / math.cos(gamma_rad)
    
    # Man's position and velocity at t1, t2, t3
    y_man_t1 = u_man_initial * t1
    y_man_t2 = y_man_t1 + u_man_initial * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    y_man_t3 = y_man_t2 # Stays still
    u_man_t3 = 0.0
    
    # Man's state at t4
    y_man_t4 = y_man_t3 + u_man_t3 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2
    u_man_t4 = u_man_t3 + a2 * (t4 - t3)
    
    # Man's position at tr
    y_man_tr = y_man_t4 + u_man_t4 * delta_t_r
    
    # The rendezvous condition is y_b(t4) = y_man_tr
    # v * y_b_t4_factor = y_man_tr
    v_bird = y_man_tr / y_b_t4_factor
    
    # Step 4: Calculate the state of man and bird at t5 = 23s
    # Man's state at t5
    delta_t_53 = t5 - t3
    y_man_t5 = y_man_t3 + 0.5 * a2 * delta_t_53**2
    u_man_t5 = a2 * delta_t_53
    
    # Bird's state at t5 (during its upward flight segment)
    x_bird_t5 = v_bird * x_b_t4_factor
    y_bird_t5 = v_bird * y_b_t4_factor
    z_bird_t5 = v_bird * (t5 - t3) # z(t3)=0 in this segment
    
    # Step 5: Analyze the final leg from t5 to t6
    delta_t_65 = t6 - t5
    
    # The bird's new velocity components to reach (0, y_final, 0)
    v_x_new = -x_bird_t5 / delta_t_65
    v_z_new = -z_bird_t5 / delta_t_65
    
    # The new northward velocity component, from v_bird^2 = vx^2+vy^2+vz^2
    v_y_new_sq = v_bird**2 - v_x_new**2 - v_z_new**2
    # The bird moves northward, so we take the positive root
    v_y_new = math.sqrt(v_y_new_sq)
    
    # Final meeting position on the y-axis
    y_final = y_bird_t5 + v_y_new * delta_t_65
    
    # Step 6: Calculate the man's final acceleration a3
    # y_final = y_man_t5 + u_man_t5*dt + 0.5*a3*dt^2
    a3 = (2 * (y_final - y_man_t5 - u_man_t5 * delta_t_65)) / (delta_t_65**2)
    
    # Output the final equation
    print(f"The man's final acceleration a3 is calculated using the kinematic equation:")
    print(f"y_final = y_man_t5 + u_man_t5 * (t6 - t5) + 0.5 * a3 * (t6 - t5)^2")
    print(f"{y_final:.2f} = {y_man_t5:.2f} + {u_man_t5:.2f} * {delta_t_65:.2f} + 0.5 * a3 * {delta_t_65:.2f}^2")
    print(f"Solving for a3 gives:")
    print(f"a3 = (2 * ({y_final:.2f} - {y_man_t5:.2f} - {u_man_t5:.2f} * {delta_t_65:.2f})) / {delta_t_65:.2f}^2")
    print(f"a3 = {a3:.2f} m/s^2")
    
    # Return the final numerical answer for the platform
    return a3

# Execute the function and print the final answer in the required format.
final_acceleration = solve_acceleration()
print(f"<<<{final_acceleration:.2f}>>>")
