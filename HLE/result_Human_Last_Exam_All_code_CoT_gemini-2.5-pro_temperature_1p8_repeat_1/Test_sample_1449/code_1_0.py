import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a free-moving wedge.
    """
    # --- Parameters from the problem ---
    m = 100 / 1000  # Convert 100 g to kg
    M = 10.0         # kg
    theta_deg = 30.0 # degrees
    h = 2.0          # m
    mu = 1.0 / 2.0     # Coefficient of friction
    g = 10.0         # m/s^2

    # --- Step 1: Preliminary calculations ---
    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    sint = math.sin(theta_rad)
    cost = math.cos(theta_rad)
    
    print("--- Physics Problem Solution ---")
    print(f"Given parameters:")
    print(f"  Block mass (m): {m} kg")
    print(f"  Wedge mass (M): {M} kg")
    print(f"  Wedge angle (theta): {theta_deg} degrees")
    print(f"  Height (h): {h} m")
    print(f"  Friction coefficient (mu): {mu}")
    print(f"  Gravity (g): {g} m/s^2\n")
    
    # --- Step 2: Solve for the acceleration of the wedge (A) ---
    # The formula for A is derived from the equations of motion for both bodies.
    # N*(sin(t) + mu*cos(t)) = M*A
    # N = m*g*cos(t) + m*A*sin(t)
    # Solving for A:
    # A = (m*g*cos(t)*(sin(t) + mu*cos(t))) / (M - m*sin(t)^2 - m*mu*sin(t)cos(t))
    numerator_A = m * g * cost * (sint + mu * cost)
    denominator_A = M - m * (sint**2) - m * mu * sint * cost
    A = numerator_A / denominator_A
    
    # --- Step 3: Solve for the relative acceleration of the block (a_rel) ---
    # The formula for a_rel is derived from the forces on the block in the wedge's frame.
    # m*a_rel = m*g*sin(t) + m*A*cos(t) - mu*N
    # a_rel = g*(sin(t) - mu*cos(t)) + A*(cos(t) - mu*sin(t))
    a_rel = g * (sint - mu * cost) + A * (cost - mu * sint)
    
    # --- Step 4: Calculate distance and time ---
    # The distance 'd' the block slides along the incline
    d = h / sint
    
    # Use kinematics: d = 0.5 * a_rel * t^2, solve for t
    if a_rel > 0:
        time = math.sqrt(2 * d / a_rel)
    else:
        print("Error: Block does not slide down (non-positive relative acceleration).")
        return

    # --- Step 5: Print results step-by-step ---
    print("--- Calculation Steps ---")
    print(f"1. Acceleration of the wedge (A): {A:.4f} m/s^2")
    print(f"2. Relative acceleration of the block (a_rel): {a_rel:.4f} m/s^2")
    print(f"3. Distance to slide along incline (d):")
    print(f"   d = h / sin(theta) = {h} / sin({theta_deg}) = {d:.4f} m")
    
    print("\n4. Final time calculation (t):")
    print(f"   Using the kinematic equation: d = (1/2) * a_rel * t^2")
    print(f"   Solving for t: t = sqrt(2 * d / a_rel)")
    print(f"   t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
    print(f"   t = {time:.4f} seconds")

    print(f"\n<<< {time:.4f} >>>")

solve_slide_time()