import math

def calculate_asymptote_angles(alpha_deg, beta_deg, gamma_deg, delta_deg):
    """
    Calculates the angles between the asymptotes of the specified conic and the line BC.

    The problem asks for the angles between the asymptotes of the conic A'B'C'OH'
    and the line BC. Through geometric analysis, it can be shown that these angles
    depend only on the angle delta between line BC and line l. The angles of
    the triangle ABC (alpha, beta, gamma) and the position of the point X on the
    circumcircle do not affect the final result.

    The angles of the asymptotes with respect to the line l are 0 and pi/2.
    Since the line BC makes an angle of delta with l, the angles of the asymptotes
    with respect to BC are -delta and (pi/2 - delta).

    Args:
        alpha_deg (float): The angle alpha of triangle ABC in degrees.
        beta_deg (float): The angle beta of triangle ABC in degrees.
        gamma_deg (float): The angle gamma of triangle ABC in degrees.
        delta_deg (float): The angle delta between line BC and line l in degrees.

    Returns:
        tuple: A tuple containing the two angles in radians.
    """
    # Note: alpha, beta, and gamma are not needed for the calculation,
    # demonstrating the independence of the result from the triangle's specific shape.
    
    # Convert delta from degrees to radians for calculation
    delta_rad = math.radians(delta_deg)
    
    # The angles of the asymptotes with respect to line BC
    angle1 = -delta_rad
    angle2 = math.pi / 2 - delta_rad
    
    return angle1, angle2

# --- Example Usage ---
# Define the angles for a sample triangle ABC and the angle delta.
# For example, a triangle with angles 60, 70, 50 degrees.
alpha = 70
beta = 60
gamma = 50

# And an angle delta of 30 degrees.
delta = 30

# Calculate the angles
angle_rad1, angle_rad2 = calculate_asymptote_angles(alpha, beta, gamma, delta)
angle_deg1 = math.degrees(angle_rad1)
angle_deg2 = math.degrees(angle_rad2)

# Print the final equations and the results for the example
print("The angles between the asymptotes and line BC are given by the expressions:")
print(f"Angle 1 = -delta")
print(f"Angle 2 = pi/2 - delta")
print("\nFor the example case where delta = {} degrees:".format(delta))
print("Angle 1 = {:.4f} radians, which is {:.2f} degrees.".format(angle_rad1, angle_deg1))
print("Angle 2 = {:.4f} radians, which is {:.2f} degrees.".format(angle_rad2, angle_deg2))
