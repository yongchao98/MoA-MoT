import math

def calculate_center_of_mass():
    """
    Calculates the horizontal and vertical coordinates of the center of mass of the candy string.
    
    Assumptions based on the problem description:
    1.  The problem can be analyzed in a 2D cross-section (e.g., the y-z plane).
        The pumpkin shape is a quarter circle of radius R. The z-axis is vertical, and the y-axis is horizontal.
    2.  The string lies on the arc from the top (y=0, z=R) to the edge (y=R, z=0).
    3.  From the edge, the string hangs vertically. The length of this hanging part is determined
        by static equilibrium.
    4.  The equilibrium analysis shows that the length of the hanging part is equal to the radius R.
    
    Let R=1 for the calculation of raw numbers.
    The total mass of the string is proportional to its total length: L_total = L_arc + L_hang = (pi*R/2) + R.
    
    The center of mass coordinates are calculated as:
    Y_cm = (integral(y*dm)_arc + integral(y*dm)_hang) / M_total
    Z_cm = (integral(z*dm)_arc + integral(z*dm)_hang) / M_total
    
    This leads to the following simplified formulas:
    Horizontal coordinate (Y_cm) = 4 * R / (pi + 2)
    Vertical coordinate (Z_cm) = R / (pi + 2)
    """
    
    # Set the radius R to 1 to get the "raw number" for the coordinates.
    R = 1.0
    
    # The value of pi from the math module.
    pi_val = math.pi
    
    # Calculate the common denominator for the coordinate expressions.
    denominator = pi_val + 2
    
    # Calculate the horizontal coordinate (along the y-axis in our model).
    # Expression: 4 * R / (pi + 2)
    horizontal_coord = 4 * R / denominator
    
    # Calculate the vertical coordinate (along the z-axis in our model).
    # Expression: R / (pi + 2)
    vertical_coord = R / denominator
    
    # Print the final coordinates separated by a comma.
    print(f"{horizontal_coord},{vertical_coord}")

calculate_center_of_mass()