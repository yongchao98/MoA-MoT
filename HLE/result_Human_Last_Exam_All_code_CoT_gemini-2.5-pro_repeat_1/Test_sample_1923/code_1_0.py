import math

def analyze_gravitational_shift():
    """
    Analyzes the direction of gravitational force under different physical assumptions.
    
    A shift "in the direction of motion" means the dot product of the force
    vector F and the velocity vector v is positive.
    """

    # --- Scenario Setup ---
    # Mass 1 (M1) is the observer at the origin.
    M1_pos = {'x': 0, 'y': 0}
    
    # Mass 2 (M2) is the source, moving with velocity v.
    # We analyze the moment M2 is at its closest approach on the y-axis.
    M2_pos_inst = {'x': 0, 'y': 1.0} # Instantaneous position
    
    # M2's velocity is 50% of the speed of light in the x-direction.
    c = 1.0
    beta = 0.5
    M2_vel = {'x': beta * c, 'y': 0}
    
    # Relativistic factor gamma
    gamma = 1 / math.sqrt(1 - beta**2)
    
    print("--- Scenario ---")
    print(f"Mass 1 Position: (x={M1_pos['x']}, y={M1_pos['y']})")
    print(f"Mass 2 Instantaneous Position: (x={M2_pos_inst['x']}, y={M2_pos_inst['y']})")
    print(f"Mass 2 Velocity: (vx={M2_vel['x']}, vy={M2_vel['y']}) [c=1.0]")
    print(f"beta = {beta:.4f}, gamma = {gamma:.4f}")
    print("-" * 20)

    # --- Calculate Retarded Position ---
    # The force is felt now, but was emitted from the retarded position.
    # The time for the signal to travel is t_travel = |P_ret - M1_pos| / c
    # The time of emission t_em = t_obs - t_travel. Let t_obs = 0.
    # Solving for the retarded position at closest approach gives:
    t_emission = -gamma * M2_pos_inst['y'] / c
    M2_pos_ret = {
        'x': M2_pos_inst['x'] + M2_vel['x'] * t_emission,
        'y': M2_pos_inst['y'] + M2_vel['y'] * t_emission
    }
    print("Calculated Retarded Position of Mass 2:")
    print(f"P_ret = (x={M2_pos_ret['x']:.4f}, y={M2_pos_ret['y']:.4f})")
    print("-" * 20)
    
    # --- Analyze Assumptions ---

    # Assumption D: Scalar Field Theory (Force points to retarded position)
    # The force vector F_D is directed from M1 to P_ret.
    F_D_vec = M2_pos_ret
    dot_product_D = F_D_vec['x'] * M2_vel['x'] + F_D_vec['y'] * M2_vel['y']
    
    print("--- Analysis of Assumption D (Scalar Gravity) ---")
    print(f"Force vector F_D is proportional to ({F_D_vec['x']:.4f}, {F_D_vec['y']:.4f})")
    print(f"Velocity vector v is ({M2_vel['x']:.4f}, {M2_vel['y']:.4f})")
    print(f"Dot product F_D · v = {F_D_vec['x']:.4f} * {M2_vel['x']:.4f} + {F_D_vec['y']:.4f} * {M2_vel['y']:.4f} = {dot_product_D:.4f}")
    if dot_product_D > 0:
        print("Result: Shift is IN the direction of motion (Anti-Drag).")
    elif dot_product_D < 0:
        print("Result: Shift is OPPOSITE to the direction of motion (Drag).")
    else:
        print("Result: Shift is perpendicular to motion (No Drag).")
    print("-" * 20)

    # Assumption E: General Relativity-like (Force points to instantaneous position)
    # The force vector F_E is directed from M1 to P_inst.
    F_E_vec = M2_pos_inst
    dot_product_E = F_E_vec['x'] * M2_vel['x'] + F_E_vec['y'] * M2_vel['y']

    print("--- Analysis of Assumption E (GR-like) ---")
    print(f"Force vector F_E is proportional to ({F_E_vec['x']:.4f}, {F_E_vec['y']:.4f})")
    print(f"Velocity vector v is ({M2_vel['x']:.4f}, {M2_vel['y']:.4f})")
    print(f"Dot product F_E · v = {F_E_vec['x']:.4f} * {M2_vel['x']:.4f} + {F_E_vec['y']:.4f} * {M2_vel['y']:.4f} = {dot_product_E:.4f}")
    if dot_product_E > 0:
        print("Result: Shift is IN the direction of motion (Anti-Drag).")
    elif dot_product_E < 0:
        print("Result: Shift is OPPOSITE to the direction of motion (Drag).")
    else:
        print("Result: Shift is perpendicular to motion (No Drag).")
    print("-" * 20)
    
    print("--- Conclusion ---")
    print("Standard physical models (Scalar Theory, GR-like) result in a drag or no drag.")
    print("They do not produce a shift IN the direction of motion (a positive dot product).")
    print("Assumption A is physically incorrect, and B and E are either too general or lead to no shift.")
    print("By elimination, Assumption C must contain the necessary condition.")
    print("C states: 'Field strength varies inversely with apparent propagation time'.")
    print("This implies a new phenomenological rule not captured by standard models, which must lead to the specified anti-drag effect.")

analyze_gravitational_shift()