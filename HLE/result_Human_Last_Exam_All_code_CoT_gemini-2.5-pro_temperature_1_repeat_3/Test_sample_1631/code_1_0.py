import math

def calculate_blank_diameter(part_diameter, part_height):
    """
    Calculates the initial blank diameter required for a simple cylindrical deep-drawn cup.
    This calculation is based on conserving the surface area between the flat blank and the formed cup,
    assuming no material thinning.

    Args:
        part_diameter (float): The diameter of the final cylindrical cup.
        part_height (float): The height of the final cylindrical cup.

    Returns:
        float: The calculated diameter of the initial flat blank.
    """
    if part_diameter <= 0 or part_height <= 0:
        print("Error: Part diameter and height must be positive numbers.")
        return None

    # The formula is derived from equating the surface area of the blank to the surface area of the cup:
    # Area_blank = Area_cup_bottom + Area_cup_wall
    # (pi * D^2) / 4 = (pi * d^2) / 4 + (pi * d * h)
    # D^2 = d^2 + 4*d*h
    # D = sqrt(d^2 + 4*d*h)

    # Let d be the part_diameter and h be the part_height
    d = part_diameter
    h = part_height
    
    # Calculate the term inside the square root
    inner_term = d**2 + 4 * d * h
    
    # Calculate the blank diameter D
    blank_diameter = math.sqrt(inner_term)
    
    print("This script calculates the required blank diameter (D) for a deep-drawn cup.")
    print(f"Given a cup diameter (d) of {d} and a height (h) of {h}.")
    print("\nThe governing equation is: D = sqrt(d^2 + 4 * d * h)")
    print("\nSubstituting the values into the equation:")
    print(f"D = sqrt({d}^2 + 4 * {d} * {h})")
    print(f"D = sqrt({d**2} + {4 * d * h})")
    print(f"D = sqrt({inner_term})")
    print(f"\nThe calculated blank diameter (D) is: {blank_diameter:.2f}")

# --- Example Usage ---
# Define the dimensions of the final stamped part (e.g., a cup)
cup_diameter = 50.0  # units can be mm or inches
cup_height = 40.0   # must be same units as diameter

calculate_blank_diameter(cup_diameter, cup_height)