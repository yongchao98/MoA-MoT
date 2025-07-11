import math

def calculate_shaded_area():
    """
    Calculates the total area of the shaded regions based on the provided image.
    """
    # Step 1: Define the given values from the diagram.
    radius = 34
    tangent_length = 17

    # Step 2: Calculate the area of the quadrilateral OATB.
    # This is twice the area of one of the right-angled triangles.
    # Area = 2 * (1/2 * radius * tangent_length) = radius * tangent_length
    area_quadrilateral = radius * tangent_length

    # Step 3: Calculate the area of the circular sector OAB.
    # First, find the half-angle 'alpha' in radians using the arctan function.
    # tan(alpha) = tangent_length / radius
    alpha = math.atan(tangent_length / radius)
    
    # The total angle of the sector is 2 * alpha.
    # The area of the sector is (1/2) * radius^2 * (2 * alpha) = radius^2 * alpha.
    area_sector = (radius**2) * alpha

    # Step 4: Calculate the total shaded area by subtracting the sector area
    # from the quadrilateral area.
    total_shaded_area = area_quadrilateral - area_sector

    # Step 5: Print the detailed calculation process.
    print("The total shaded area is the area of the quadrilateral minus the area of the circular sector.")
    print("\n--- Calculation Breakdown ---")
    
    # Print the equation for the quadrilateral's area
    print("\n1. Area of the Quadrilateral:")
    print(f"   Equation: radius * tangent_length")
    print(f"   Calculation: {radius} * {tangent_length} = {area_quadrilateral}")

    # Print the equation for the sector's area
    print("\n2. Area of the Circular Sector:")
    print(f"   Equation: radius^2 * arctan(tangent_length / radius)")
    print(f"   Calculation: {radius}^2 * arctan({tangent_length} / {radius})")
    print(f"   = {radius**2} * arctan(0.5)")
    print(f"   = {area_sector}")

    # Print the final equation for the total shaded area
    print("\n3. Total Shaded Area:")
    print(f"   Equation: Area of Quadrilateral - Area of Sector")
    print(f"   Calculation: {area_quadrilateral} - {area_sector}")
    print(f"   Final Result: {total_shaded_area}")

# Run the calculation
calculate_shaded_area()