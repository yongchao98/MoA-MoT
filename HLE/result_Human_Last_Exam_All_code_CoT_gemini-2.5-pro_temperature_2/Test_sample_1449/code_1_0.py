import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 100 / 1000  # mass of block in kg
    M = 10.0        # mass of wedge in kg
    theta_deg = 30.0# angle in degrees
    h = 2.0         # height in m
    mu = 1.0/2.0    # coefficient of friction
    g = 10.0        # acceleration due to gravity in m/s^2

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sint = math.sin(theta_rad)
    cost = math.cos(theta_rad)
    
    # --- Step 1: Calculate the acceleration of the wedge (A_x) ---
    # Derived from Newton's laws on both bodies.
    # A_x = (m*g*cos(theta)*(sin(theta) + mu*cos(theta))) / (M - m*sin(theta)^2 - m*mu*sin(theta)*cos(theta))
    
    ax_numerator = m * g * cost * (sint + mu * cost)
    ax_denominator = M - m * sint**2 - m * mu * sint * cost
    
    A_x = ax_numerator / ax_denominator
    
    # --- Step 2: Calculate the block's acceleration relative to the wedge (a_rel) ---
    # Derived from forces on the block in the non-inertial frame.
    # a_rel = g*(sin(theta) - mu*cos(theta)) - A_x*(mu*sin(theta) + cos(theta))

    term1 = g * (sint - mu * cost)
    term2 = A_x * (mu * sint + cost)
    a_rel = term1 - term2

    # --- Step 3: Calculate the distance (L) and time (t) ---
    # Distance L along the incline
    L = h / sint
    
    # Time t to travel distance L from rest with acceleration a_rel
    # L = 0.5 * a_rel * t^2 => t = sqrt(2 * L / a_rel)
    if a_rel > 0:
        time = math.sqrt(2 * L / a_rel)
    else:
        # If acceleration is not positive, the block won't slide down
        time = float('inf')

    # --- Step 4: Output the results clearly ---
    print(f"Given parameters:")
    print(f"Block mass m = {m} kg")
    print(f"Wedge mass M = {M} kg")
    print(f"Wedge angle theta = {theta_deg} degrees")
    print(f"Height h = {h} m")
    print(f"Friction coefficient mu = {mu}")
    print(f"Gravity g = {g} m/s^2")
    print("-" * 30)

    print("Step 1: The acceleration of the wedge is calculated.")
    print(f"A_x = {A_x:.4f} m/s^2")
    print("-" * 30)

    print("Step 2: The acceleration of the block relative to the wedge is calculated.")
    print(f"a_rel = g*(sin({theta_deg}) - {mu}*cos({theta_deg})) - A_x*(mu*sin({theta_deg}) + cos({theta_deg}))")
    print(f"a_rel = {term1:.4f} - {A_x:.4f} * {(mu * sint + cost):.4f}")
    print(f"a_rel = {a_rel:.4f} m/s^2")
    print("-" * 30)
    
    print("Step 3: The time is found using kinematics.")
    print(f"The distance along the incline is L = h / sin({theta_deg}) = {h} / {sint:.3f} = {L:.4f} m.")
    print("The final equation for time is t = sqrt(2 * L / a_rel).")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print("-" * 30)
    print(f"Final Answer: The time for the block to reach the bottom is {time:.4f} seconds.")

    return time

if __name__ == '__main__':
    final_time = solve_sliding_time()
    # Final answer in the requested format
    # The printed output already explains everything, here is just the required marker.
    # The result is rounded to 4 decimal places in the printout, so we will use that for the final answer.
    # Using 3.7163 from the actual calculation.
    print(f'<<<{final_time:.4f}>>>')
