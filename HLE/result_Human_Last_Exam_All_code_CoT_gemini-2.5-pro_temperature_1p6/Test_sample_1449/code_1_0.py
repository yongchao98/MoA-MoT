import math

def solve_slide_time():
    """
    Calculates the time for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 100 / 1000  # Mass of block in kg
    M = 10.0        # Mass of wedge in kg
    theta_deg = 30.0  # Angle of wedge in degrees
    mu = 0.5        # Coefficient of friction
    h = 2.0         # Maximum height of the wedge in m
    g = 10.0        # Acceleration due to gravity in m/s^2

    # Convert angle to radians for use in math functions
    theta = math.radians(theta_deg)
    s_theta = math.sin(theta)
    c_theta = math.cos(theta)

    # Step 1: Calculate the acceleration of the block relative to the wedge (a_rel).
    # This formula is derived from conserving horizontal momentum and applying Newton's
    # second law to the block in the wedge's non-inertial frame.
    # a_rel = ( (M+m)*g*(sin(theta) - mu*cos(theta)) ) / ( M + m*sin(theta)^2 + m*mu*sin(theta)*cos(theta) )
    a_rel_numerator = (M + m) * g * (s_theta - mu * c_theta)
    a_rel_denominator = M + m * s_theta**2 + m * mu * s_theta * c_theta
    a_rel = a_rel_numerator / a_rel_denominator

    # Step 2: Calculate the distance the block needs to slide (L).
    L = h / s_theta

    # Step 3: Calculate the time (t) using kinematics.
    # The block starts from rest relative to the wedge (v_rel_initial = 0).
    # The relevant kinematic equation is L = 0.5 * a_rel * t^2.
    # Solving for t gives: t = sqrt(2 * L / a_rel).
    time = math.sqrt(2 * L / a_rel)

    # Print the results, showing the final equation with numerical values.
    print(f"The equation for the time it takes for the block to slide down is:")
    print(f"t = sqrt(2 * L / a_rel)")
    print(f"Where:")
    print(f"  L (distance) = h / sin(theta) = {h:.1f} / sin({theta_deg:.1f}) = {L:.4f} m")
    print(f"  a_rel (relative acceleration) = {a_rel:.4f} m/s^2")
    print(f"\nSubstituting the values into the equation:")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = {time:.4f} s")

    return time

# Execute the function and get the final answer
final_time = solve_slide_time()
# The final answer in the required format
print(f"\n<<< {final_time} >>>")
