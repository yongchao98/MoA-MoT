import math

def solve():
    """
    Solves the physics problem to find the man's final acceleration.
    """
    # --- Given constants ---
    u = 1.5  # m/s
    t0 = 0
    t1 = 4   # s
    a1 = -0.15 # m/s^2
    t3 = 15  # s
    a2 = 0.25  # m/s^2
    gamma_deg = 40
    t5 = 23  # s
    t6 = 40  # s

    # --- Helper functions for trigonometry in degrees ---
    def sind(deg): return math.sin(math.radians(deg))
    def cosd(deg): return math.cos(math.radians(deg))

    # --- Step 1: Man's timeline and positions ---
    y_m_t1 = u * (t1 - t0) # Position at t1=4s
    
    t_decel = u / abs(a1) # Time to decelerate from u to 0
    t2 = t1 + t_decel     # Time when man stops
    
    y_m_t2 = y_m_t1 + u * t_decel + 0.5 * a1 * t_decel**2 # Position at t2
    y_m_t3 = y_m_t2 # Man is still from t2 to t3=15s
    
    t_accel = u / a2 # Time to accelerate from 0 to u
    t4 = t3 + t_accel  # Time when man reaches speed u again
    
    y_m_t4 = y_m_t3 + 0.5 * a2 * t_accel**2 # Position at t4
    
    # --- Step 2, 3, 4: Solve for bird's speed v ---
    # From the rendezvous conditions, a quadratic equation for r = v_xy/v is derived:
    # 13.268 * r^2 + 7.218 * r - 6.032 = 0
    A, B, C = 13.268, 7.218, -6.032
    r = (-B + math.sqrt(B**2 - 4 * A * C)) / (2 * A)
    
    # With r, solve for v using the second rendezvous condition
    c_gamma = cosd(gamma_deg)
    s_gamma = sind(gamma_deg)
    v_numerator = y_m_t4 + u * (4 * r + 10 / c_gamma)
    v_denominator = (4 * r * s_gamma + 1)
    v = v_numerator / v_denominator

    # --- Step 5: State at t5 = 23s ---
    # Man's position at t5
    y_m_t5 = y_m_t4 + u * (t5 - t4)
    v_m_t5 = u

    # Bird's position at t5
    v_xy = r * v
    v_x1 = v_xy * c_gamma
    v_y1 = v_xy * s_gamma
    v_z1 = math.sqrt(v**2 - v_xy**2)
    
    dt_seg1 = t1-t0
    dt_seg2 = t2-t1
    dt_seg3 = t3-t2
    dt_seg4 = t4-t3
    
    x_b_t4 = dt_seg1 * v_x1 + dt_seg2 * v
    y_b_t4 = dt_seg1 * v_y1 + dt_seg3 * v
    z_b_t4 = dt_seg1 * v_z1 + dt_seg4 * v

    v_bx5_planned = -v * c_gamma
    v_bz5_planned = -v * s_gamma

    dt_to_t5 = t5 - t4
    x_b_t5 = x_b_t4 + v_bx5_planned * dt_to_t5
    y_b_t5 = y_b_t4 # No y-motion in this segment
    z_b_t5 = z_b_t4 + v_bz5_planned * dt_to_t5
    
    # --- Step 6 & 7: Solve for a3 ---
    dt_final = t6 - t5
    
    # Bird displacement components from t5 to t6
    dx_bird = 0 - x_b_t5
    dz_bird = 0 - z_b_t5
    
    # Squared distance bird must travel in x and z
    dist_sq_xz_bird = dx_bird**2 + dz_bird**2
    
    # Total squared distance bird travels
    total_dist_sq_bird = (v * dt_final)**2
    
    # This implies the squared displacement in y is
    dy_bird_sq = total_dist_sq_bird - dist_sq_xz_bird
    # Bird moves northward, so displacement dy is positive
    dy_bird = math.sqrt(dy_bird_sq)

    # Man's displacement in y
    # dy_man = y_meet - y_m_t5 = v_m_t5 * dt_final + 0.5 * a3 * dt_final**2
    # The y displacements must be equal, dy_man == dy_bird
    # dy_bird = v_m_t5 * dt_final + 0.5 * a3 * dt_final**2
    
    a3_numerator = dy_bird - (v_m_t5 * dt_final)
    a3_denominator = 0.5 * dt_final**2
    
    a3 = a3_numerator / a3_denominator
    
    print("The final leg of the journey begins at t_5 = 23 s and ends at t_6 = 40 s.")
    print("The final acceleration, a_3, is found by equating the man's and the bird's displacement in the North (y) direction during this time.")
    print("\nFinal equation for the acceleration a_3:")
    print(f"a_3 = ( Bird_y_displacement - (Man_initial_speed * time) ) / ( 0.5 * time^2 )")
    print("\nWhere:")
    print(f"Bird_y_displacement = {dy_bird:.4f} m")
    print(f"Man_initial_speed = {v_m_t5:.4f} m/s")
    print(f"time = {dt_final:.4f} s")
    print("\nPlugging in the numbers:")
    print(f"a_3 = ( {dy_bird:.4f} - ({v_m_t5:.4f} * {dt_final:.4f}) ) / ( 0.5 * {dt_final:.4f}^2 )")
    print(f"a_3 = ( {dy_bird:.4f} - {v_m_t5 * dt_final:.4f} ) / {a3_denominator:.4f}")
    print(f"a_3 = {a3_numerator:.4f} / {a3_denominator:.4f}")
    print(f"a_3 = {a3:.4f} m/s^2")
    
    # The final answer in the required format.
    print(f"\n<<<{a3:.4f}>>>")

solve()