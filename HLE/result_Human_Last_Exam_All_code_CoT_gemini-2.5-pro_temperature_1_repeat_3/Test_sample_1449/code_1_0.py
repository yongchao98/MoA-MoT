import math

def solve_sliding_time():
    """
    Calculates the time for a block to slide down a moving wedge with friction.
    """
    # Given parameters
    m = 100 / 1000  # Convert grams to kilograms
    M = 10.0        # kg
    theta_deg = 30.0
    h = 2.0         # m
    mu = 1.0 / 2.0
    g = 10.0        # m/s^2

    # The initial velocity v0 is not needed to calculate the time for the block
    # to slide down, as the accelerations are independent of the initial velocity.

    print("--- Step 1: Define constants and calculate sliding distance L ---")
    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    s_theta = math.sin(theta_rad)
    c_theta = math.cos(theta_rad)

    # Calculate the distance L the block slides along the wedge's surface
    L = h / s_theta
    print(f"Wedge angle: {theta_deg} degrees")
    print(f"Wedge height h: {h:.4f} m")
    print(f"Sliding distance L = h / sin(theta) = {L:.4f} m\n")

    print("--- Step 2: Calculate the acceleration of the wedge (A) and the block (a_rel) ---")
    # The equations of motion for the system are solved to find the accelerations.
    # The derivation leads to the following formulas for A and a_rel.

    # Formula for the wedge's horizontal acceleration, A
    A_numerator = m * g * c_theta * (s_theta + mu * c_theta)
    A_denominator = M - m * s_theta * (s_theta + mu * c_theta)
    A = A_numerator / A_denominator
    print(f"Calculated acceleration of the wedge A: {A:.4f} m/s^2")

    # Formula for the block's acceleration relative to the wedge, a_rel
    a_rel = g * (s_theta - mu * c_theta) - A * (mu * s_theta + c_theta)
    print(f"Calculated relative acceleration of the block a_rel: {a_rel:.4f} m/s^2\n")

    print("--- Step 3: Calculate the time t ---")
    # Check if the block will slide down at all
    if a_rel <= 0:
        print("The block does not slide down the wedge (acceleration is not positive).")
        t = float('inf')
    else:
        # Use the kinematic equation: L = (1/2) * a_rel * t^2
        t = math.sqrt(2 * L / a_rel)
        
        # Output the final equation with the calculated numbers
        print("The final equation for time t is: t = sqrt( (2 * L) / a_rel )")
        print(f"Plugging in the numbers: t = sqrt( (2 * {L:.4f}) / {a_rel:.4f} )")
        print(f"\nThe final time is {t:.4f} seconds.")

    return t

# Run the solver and get the final answer
final_time = solve_sliding_time()

# Output the final answer in the required format
print(f"\n<<< {final_time:.4f} >>>")
