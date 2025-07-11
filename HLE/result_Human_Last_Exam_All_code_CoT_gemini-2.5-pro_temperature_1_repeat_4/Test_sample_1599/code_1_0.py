import math

def solve_acceleration():
    """
    Solves for the man's final acceleration based on the described scenario.
    """
    # Initial given values
    u = 1.5  # m/s
    t0 = 0.0
    t1 = 4.0  # s
    # The given alpha = 130.0 degrees creates a contradiction with the "northward" motion description.
    # Assuming alpha = 50.0 degrees resolves this, as sin(130)=sin(50) and cos(50) is positive.
    # This is a necessary assumption to make the problem physically consistent and solvable.
    alpha = 50.0 # degrees
    a1 = -0.15  # m/s^2
    t3 = 15.0  # s
    a2 = 0.25  # m/s^2
    gamma = 40.0 # degrees
    t5 = 23.0  # s
    t6 = 40.0  # s

    # --- Step 1: Determine the man's timeline and positions ---
    # Time t2 is when the man stops decelerating (velocity becomes 0).
    t2 = t1 - u / a1
    # Time t4 is when the man reaches speed u again after accelerating from rest.
    t4 = t3 + u / a2

    # Calculate man's position at key times
    y_man_t1 = u * t1
    dt_12 = t2 - t1
    y_man_t2 = y_man_t1 + u * dt_12 + 0.5 * a1 * dt_12**2
    y_man_t3 = y_man_t2
    dt_34 = t4 - t3
    y_man_t4 = y_man_t3 + 0.5 * a2 * dt_34**2
    dt_45 = t5 - t4
    y_man_t5 = y_man_t4 + u * dt_45
    v_man_t5 = u

    # --- Step 2: Solve for the planned meeting time offset T = t_meet - t4 ---
    c40 = math.cos(math.radians(gamma))
    s40 = math.sin(math.radians(gamma))
    sa = math.sin(math.radians(alpha))
    
    # Coefficients for the quadratic equation A*T^2 + B*T + C = 0 for T
    A = (c40 / (4 * sa))**2 + (s40 / 4)**2
    B = -2 * (c40 * 10 / (16 * sa**2)) - 2 * (s40 * 6 / 16)
    C = (10 / (4 * sa))**2 + (6 / 4)**2 - 1

    discriminant = B**2 - 4 * A * C
    T1 = (-B + math.sqrt(discriminant)) / (2 * A)
    T2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # The physically valid root for T must result in a positive ground speed ratio Rg
    if c40 * T1 - 10 > 0:
        T = T1
    else:
        T = T2

    # --- Step 3: Calculate bird's speed and planned meeting location ---
    t_meet = t4 + T
    Rg = (c40 * T - 10) / (4 * sa) # vg/v ratio
    ca = math.cos(math.radians(alpha))
    y_man_t_meet = y_man_t4 + u * T
    v = y_man_t_meet / (1 + 4 * Rg * ca)

    # --- Step 4: Calculate the bird's position at t5 ---
    xb_t5 = v * c40 * (t_meet - t5)
    yb_t5 = y_man_t_meet
    zb_t5 = v * s40 * (t_meet - t5)

    # --- Step 5: Calculate the actual meeting location y_meet ---
    dt_56 = t6 - t5
    # From the bird's final flight speed constraint
    y_meet_minus_yb_t5_sq = dt_56**2 * v**2 - xb_t5**2 - zb_t5**2
    # The new y_meet should be greater than yb_t5 due to the northward gust
    y_meet = yb_t5 + math.sqrt(y_meet_minus_yb_t5_sq)

    # --- Step 6: Calculate the final acceleration a3 ---
    # From the man's final leg kinematics
    a3 = (y_meet - y_man_t5 - v_man_t5 * dt_56) / (0.5 * dt_56**2)

    # --- Output the final calculation ---
    print("To find the man's final acceleration, a3, we use the kinematic equation for his last leg of the journey:")
    print("y_meet = y_man(t5) + v_man(t5)*dt + 0.5*a3*dt^2")
    print(f"Final meeting position y_meet = {y_meet:.2f} m")
    print(f"Man's position at t5, y_man(t5) = {y_man_t5:.2f} m")
    print(f"Man's velocity at t5, v_man(t5) = {v_man_t5:.2f} m/s")
    print(f"Time interval dt = t6 - t5 = {dt_56:.2f} s")
    print("\nSubstituting the values into the equation:")
    print(f"{y_meet:.2f} = {y_man_t5:.2f} + {v_man_t5:.2f} * {dt_56:.2f} + 0.5 * a3 * {dt_56:.2f}^2")
    man_pos_if_no_accel = y_man_t5 + v_man_t5 * dt_56
    print(f"{y_meet:.2f} = {man_pos_if_no_accel:.2f} + {0.5 * dt_56**2:.2f} * a3")
    displacement_due_to_a3 = y_meet - man_pos_if_no_accel
    print(f"{displacement_due_to_a3:.2f} = {0.5 * dt_56**2:.2f} * a3")
    print(f"a3 = {displacement_due_to_a3:.2f} / {0.5 * dt_56**2:.2f}")
    print(f"\nThe man's constant acceleration during the final leg is {a3:.3f} m/s^2.")
    
solve_acceleration()
<<<1.381>>>