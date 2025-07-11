import math

def calculate_center_of_mass(R):
    """
    Calculates the center of mass of a candy string on a spherical pumpkin.

    Args:
      R (float): The radius of the spherical pumpkin.

    Returns:
      tuple: A tuple containing the (x, y) coordinates of the center of mass.
    """

    # The string is composed of two parts:
    # Part 1: A quarter-circular arc of radius R.
    # Part 2: A vertical line segment of length R.
    # The linear density rho is constant and cancels out in the final calculation.

    # -- Properties of Part 1 (Quarter-circular arc) --
    # Length of the arc is 1/4 of a circle's circumference.
    # The arc goes from (0, R) to (R, 0).
    m1_factor = math.pi / 2  # Mass is proportional to rho * R * (pi/2)
    # The center of mass of a quarter-circular arc in the first quadrant is (2R/pi, 2R/pi).
    x_cm_1 = (2 * R) / math.pi
    y_cm_1 = (2 * R) / math.pi

    # -- Properties of Part 2 (Vertical line segment) --
    # Length of the hanging segment is R.
    # It hangs from (R, 0) to (R, -R).
    m2_factor = 1.0  # Mass is proportional to rho * R * 1
    # The center of mass is the midpoint of the segment.
    x_cm_2 = R
    y_cm_2 = -R / 2

    # -- Total Center of Mass Calculation --
    # Total mass is proportional to the sum of the mass factors.
    total_m_factor = m1_factor + m2_factor

    # Horizontal coordinate (X_cm)
    # The equation is: X_cm = (m1*x1 + m2*x2) / (m1 + m2)
    # Plugging in the formulas gives: X_cm = (R^2 + R^2) / (R*(pi/2 + 1)) = 4*R / (pi + 2)
    # Here are the numbers in the final simplified equation for X_cm:
    numerator_x_factor = 4
    denominator_pi_factor = 1
    denominator_const_factor = 2
    
    x_cm = (numerator_x_factor * R) / (denominator_pi_factor * math.pi + denominator_const_factor)

    # Vertical coordinate (Y_cm)
    # The equation is: Y_cm = (m1*y1 + m2*y2) / (m1 + m2)
    # Plugging in the formulas gives: Y_cm = (R^2 - R^2/2) / (R*(pi/2 + 1)) = R / (pi + 2)
    # Here are the numbers in the final simplified equation for Y_cm:
    numerator_y_factor = 1
    
    y_cm = (numerator_y_factor * R) / (denominator_pi_factor * math.pi + denominator_const_factor)

    return x_cm, y_cm

# You can set the radius R to any value.
# We will use a default value of R=1 for the calculation.
R = 1.0

# Calculate the coordinates
x_coordinate, y_coordinate = calculate_center_of_mass(R)

# Print the final raw numbers as requested, separated by a comma.
print(f"{x_coordinate},{y_coordinate}")
