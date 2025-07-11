import math

def solve_area():
    """
    Calculates the total area of the shaded regions based on the provided image.
    """
    # Define the given values
    radius = 34
    length = 17

    # The total shaded area is the sum of two identical regions.
    # The area of one region is the area of a right-angled triangle minus the area of a circular sector.

    # 1. Calculate the combined area of the two construction triangles.
    # Each triangle has legs of 'radius' and 'length'.
    # Area = 2 * (1/2 * base * height)
    area_of_two_triangles = 2 * (0.5 * radius * length)

    # 2. Calculate the combined area of the two sectors being subtracted.
    # The angle of each sector is arctan(length / radius).
    # Area = 2 * (1/2 * radius^2 * angle)
    angle_rad = math.atan(length / radius)
    area_of_two_sectors = 2 * (0.5 * (radius**2) * angle_rad)

    # 3. The total shaded area is the difference.
    total_area = area_of_two_triangles - area_of_two_sectors
    
    # Print the equation with all the numbers
    print("The total shaded area is calculated by combining two identical shaded regions.")
    print("Area of one region = Area of a right triangle - Area of a circular sector.")
    print(f"Total Area = 2 * ( (1/2 * {radius} * {length}) - (1/2 * {radius}**2 * arctan({length}/{radius})) )")
    print(f"Total Area = ({radius} * {length}) - ({radius}**2 * arctan(0.5))")
    print(f"Total Area = {int(radius * length)} - {radius**2} * arctan(0.5)")
    print(f"Total Area = {area_of_two_triangles} - {area_of_two_sectors:.4f}")
    print(f"Total Area ≈ {total_area:.4f}")

solve_area()