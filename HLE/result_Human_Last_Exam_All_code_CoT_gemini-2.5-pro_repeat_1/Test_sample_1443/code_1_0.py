import math

def get_asymptote_angles(alpha_deg, beta_deg, gamma_deg, delta_deg):
    """
    Calculates the angles between the asymptotes of the specified conic and the line BC.

    Args:
      alpha_deg (float): Angle A of triangle ABC in degrees.
      beta_deg (float): Angle B of triangle ABC in degrees.
      gamma_deg (float): Angle C of triangle ABC in degrees.
      delta_deg (float): Angle between line l and line BC in degrees.
    
    Returns:
      tuple: A tuple containing the two angles in degrees.
    """
    
    # According to the geometric properties of this configuration (IMO SL 2011 G8),
    # the conic is a rectangular hyperbola. One of its axes is parallel to the line l.
    # The asymptotes of a rectangular hyperbola are at ±45 degrees (or pi/4 radians)
    # to its axes.
    # Therefore, the asymptotes make angles of delta ± 45 degrees with the line BC.

    # Convert delta from degrees to radians for calculation
    delta_rad = math.radians(delta_deg)
    
    # Angle difference is pi/4 radians
    angle_diff_rad = math.pi / 4
    
    # Calculate the two asymptote angles in radians relative to BC
    angle1_rad = delta_rad + angle_diff_rad
    angle2_rad = delta_rad - angle_diff_rad
    
    # Convert back to degrees for the output
    angle1_deg = math.degrees(angle1_rad)
    angle2_deg = math.degrees(angle2_rad)
    
    return angle1_deg, angle2_deg

# Example values for the angles, not used in the final formula but required by the function signature.
# The sum of triangle angles should be 180.
alpha = 60
beta = 70
gamma = 50

# The angle between line l and BC
delta = 30

# Calculate the angles
angle1, angle2 = get_asymptote_angles(alpha, beta, gamma, delta)

# Print the final result in a descriptive way
# The formula for the angles is delta +/- 45 degrees.
# We will show the calculation explicitly.
print(f"Given the angles of the triangle as alpha={alpha}, beta={beta}, gamma={gamma} degrees,")
print(f"and the angle between line l and BC as delta={delta} degrees.")
print("The problem describes a rectangular hyperbola whose axes are rotated by the angle delta with respect to BC.")
print("The asymptotes of a rectangular hyperbola are at +/- 45 degrees from its axes.")
print("Therefore, the angles of the asymptotes with respect to the line BC are:")
print(f"Angle 1 = delta + 45 = {delta} + 45 = {angle1:.2f} degrees")
print(f"Angle 2 = delta - 45 = {delta} - 45 = {angle2:.2f} degrees")
print("\nThe final expression for the angles in terms of delta is delta +/- 45 degrees (or delta +/- pi/4 radians).")
print("The angles do not depend on alpha, beta, or gamma.")
