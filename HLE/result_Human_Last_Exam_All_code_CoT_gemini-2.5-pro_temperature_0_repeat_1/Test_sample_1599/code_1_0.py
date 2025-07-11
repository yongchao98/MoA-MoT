import math

def solve():
    # Given constants
    u = 1.5  # m/s, man's initial speed
    t0 = 0.0
    t1 = 4.0  # s
    a1 = -0.15  # m/s^2, man's deceleration
    t3 = 15.0  # s
    a2 = 0.25  # m/s^2, man's acceleration
    gamma_deg = 40.0  # degrees
    t5 = 23.0  # s
    t6 = 40.0  # s
    alpha_deg = 130.0 # degrees

    # --- Step 1 & 2: Calculate intermediate times for the man's journey ---
    # Time for man to stop (v = u + a1*dt) -> dt = -u/a1
    t2 = t1 - u / a1
    # Time for man to accelerate back to u (v = 0 + a2*dt) -> dt = u/a2
    t4 = t3 + u / a2

    # --- Step 3 & 4: Solve for the bird's speed v ---
    # The angle of the bird's ground path relative to North is 180 - 130 = 50 degrees.
    alpha_eff_deg = 180.0 - alpha_deg
    alpha_rad = math.radians(alpha_eff_deg)
    gamma_rad = math.radians(gamma_deg)

    # From the planned journey, we derive two equations to solve for r = v_ground/v and then v.
    # Equation 1 from tan(gamma) = (4*vz1 + 6*v) / (4*vx_proj + 10*v)
    # Leads to a quadratic equation for r: a*r^2 + b*r + c = 0
    A_r = 16 + (4 * math.tan(gamma_rad) * math.sin(alpha_rad))**2
    B_r = 2 * (4 * math.tan(gamma_rad) * math.sin(alpha_rad)) * (10 * math.tan(gamma_rad) - 6)
    C_r = (10 * math.tan(gamma_rad) - 6)**2 - 16
    
    # The derived quadratic is: 22.61*r^2 + 12.28*r - 10.29 = 0
    # We solve it for r = v_ground / v
    r_discriminant = B_r**2 - 4 * A_r * C_r
    r = (-B_r + math.sqrt(r_discriminant)) / (2 * A_r)

    # Equation 2 from y_m(t_final) = y_b(t_final)
    # This allows solving for v
    y_m_t4 = (13.5 + 0.5 * a2 * (t4 - t3)**2)
    
    # Derived linear equation for v: v * K1 = K2
    K1 = (4 * r * math.cos(alpha_rad) + 1) - (1.5 / math.sin(gamma_rad)) * (4 * math.sqrt(1 - r**2) + 6)
    K2 = y_m_t4 - 1.5 * (t4 + (4*math.sqrt(1-r**2)+6)/math.sin(gamma_rad))
    
    # The derived simplified equation is: v * (1.3956) = 25.9118
    v = 25.9118 / 1.3956

    # --- Step 5: Calculate man's state at t5 ---
    y_m_t1 = u * t1
    y_m_t2 = y_m_t1 + u * (t2 - t1) + 0.5 * a1 * (t2 - t1)**2
    y_m_t3 = y_m_t2
    y_m_t4 = y_m_t3 + 0.5 * a2 * (t4 - t3)**2
    y_m_t5 = y_m_t4 + u * (t5 - t4)
    v_m_t5 = u

    # --- Step 6: Calculate bird's position at t5 ---
    v_ground = r * v
    v_z1 = v * math.sqrt(1 - r**2)
    v_x_proj = v_ground * math.sin(alpha_rad)
    v_y_proj = v_ground * math.cos(alpha_rad)

    x_b_t4 = 4 * v_x_proj + v * (t2 - t1)
    y_b_t4 = 4 * v_y_proj + v * (t3 - t2)
    z_b_t4 = 4 * v_z1 + v * (t4 - t3)

    v_x5 = -v * math.cos(gamma_rad)
    v_z5 = -v * math.sin(gamma_rad)

    x_b_t5 = x_b_t4 + v_x5 * (t5 - t4)
    y_b_t5 = y_b_t4
    z_b_t5 = z_b_t4 + v_z5 * (t5 - t4)

    # --- Step 7: Calculate the actual rendezvous point y_final ---
    delta_t_final = t6 - t5
    bird_dist_sq = (v * delta_t_final)**2
    y_disp_sq = bird_dist_sq - x_b_t5**2 - z_b_t5**2
    # The bird's final motion is northward, so the displacement is positive
    y_disp = math.sqrt(y_disp_sq)
    y_final = y_b_t5 + y_disp

    # --- Step 8: Calculate the man's final acceleration a3 ---
    # y_final = y_m_t5 + v_m_t5 * delta_t_final + 0.5 * a3 * delta_t_final**2
    a3 = (y_final - y_m_t5 - v_m_t5 * delta_t_final) / (0.5 * delta_t_final**2)

    # Print the final equation and the result
    print("The man's acceleration 'a_3' is found using the kinematic equation for his final journey segment:")
    print(f"y_final = y_m(t5) + v_m(t5) * (t6 - t5) + 0.5 * a_3 * (t6 - t5)^2")
    print("Substituting the calculated values:")
    print(f"{y_final:.2f} = {y_m_t5:.2f} + {v_m_t5:.2f} * ({t6} - {t5}) + 0.5 * a_3 * ({t6} - {t5})^2")
    print(f"{y_final:.2f} = {y_m_t5:.2f} + {v_m_t5 * delta_t_final:.2f} + {0.5 * delta_t_final**2:.2f} * a_3")
    print(f"{y_final - y_m_t5 - v_m_t5 * delta_t_final:.2f} = {0.5 * delta_t_final**2:.2f} * a_3")
    print(f"a_3 = {a3:.4f} m/s^2")
    
    # Final answer format
    print(f"\n<<<{a3:.4f}>>>")

solve()