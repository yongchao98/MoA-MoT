import math

def solve():
    """
    Calculates the time for a block to slide down a movable wedge with friction.
    """
    # Step 1: Define the given parameters
    m = 100 / 1000  # mass of the block in kg (100 g)
    M = 10.0        # mass of the wedge in kg
    theta_deg = 30.0  # angle of the wedge in degrees
    h = 2.0         # maximum height of the wedge in meters
    mu = 0.5        # coefficient of friction
    g = 10.0        # acceleration due to gravity in m/s^2

    # Step 2: Set up for calculations
    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)
    s = math.sin(theta_rad)
    c = math.cos(theta_rad)

    # Step 3: Solve for the accelerations
    # Based on the laws of motion, we derive the formulas for the acceleration
    # of the wedge (A) and the relative acceleration of the block (a_rel).

    # The acceleration of the wedge (A) can be found using:
    # A = (m*g*c*(s + mu*c)) / (M - m*s*(s + mu*c))
    A_numerator = m * g * c * (s + mu * c)
    A_denominator = M - m * s * (s + mu * c)
    A = A_numerator / A_denominator

    # The acceleration of the block relative to the wedge (a_rel) is:
    # a_rel = g*(s - mu*c) + A*(c - mu*s)
    a_rel = g * (s - mu * c) + A * (c - mu * s)

    # Step 4: Calculate the time
    # The distance 'd' the block needs to travel down the incline
    d = h / s
    
    # Using the kinematic equation d = 0.5 * a_rel * t^2
    # We solve for t = sqrt(2 * d / a_rel)
    t = math.sqrt(2 * d / a_rel)

    # Step 5: Output the results clearly
    print("The final equation for the time 't' is t = sqrt(2 * d / a_rel), where 'd' is the distance along the incline and 'a_rel' is the block's acceleration relative to the wedge.")
    print("\n--- Intermediate Calculations ---")
    print(f"Distance to travel, d = h / sin(theta) = {h:.1f} / {s:.4f} = {d:.4f} m")
    print(f"Acceleration of the wedge, A = {A:.4f} m/s^2")
    print(f"Acceleration of the block relative to the wedge, a_rel = {a_rel:.4f} m/s^2")
    
    print("\n--- Final Calculation ---")
    print(f"t = sqrt(2 * d / a_rel)")
    print(f"t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
    print(f"t = {t:.4f} s")
    
    # Print the final answer in the required format
    print(f"\n<<<The exact amount of time it takes for the block to slide all the way to the bottom of the wedge is {t} seconds.>>>")
    

solve()
<<<3.333590425713437>>>