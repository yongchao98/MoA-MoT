import math

def solve_cube_locus_problem():
    """
    This function calculates the requested percentage based on the geometry of the problem.
    The final result is independent of the cube's side length 'r'.
    """
    # Step 1: Define the components of the total length of the locus C.
    # The locus is composed of a number of identical circular arcs.
    # Number of arcs: From geometric analysis, we find there are 6 arcs.
    num_arcs = 6

    # Angle of each arc (in radians): Each arc subtends an angle of 60 degrees (pi/3 radians).
    # We can represent the angle as a fraction of pi.
    angle_numerator = 1
    angle_denominator = 3
    
    # The total length of C is num_arcs * r * (angle_numerator * pi / angle_denominator)
    # Total length = 6 * r * (pi/3) = 2 * pi * r

    # Step 2: Define the components of the normalizing factor from the problem statement.
    # We need to divide the total length by (2 * pi * r).
    divisor_coefficient = 2
    
    # Step 3: Calculate the final ratio and percentage.
    # Ratio = (6 * r * pi / 3) / (2 * pi * r)
    # The terms 'r' and 'pi' cancel out.
    ratio = (num_arcs / angle_denominator) / divisor_coefficient
    
    # Convert the ratio to a whole number percentage.
    percentage_multiplier = 100
    final_percentage = int(round(ratio * percentage_multiplier))
    
    # Step 4: Output the explanation and the final equation with its numbers.
    print("The total length of the locus C is given by the formula:")
    print("Total Length = (Number of Arcs) * r * (Arc Angle)")
    print(f"             = {num_arcs} * r * (pi / {angle_denominator})")
    print("\nThe problem asks to divide this by (2 * pi * r) and express it as a whole number percentage.")
    print("\nThe final equation with its numerical components is:")
    print(f"Percentage = (({num_arcs} * r * pi) / {angle_denominator}) / ({divisor_coefficient} * pi * r) * {percentage_multiplier}")
    print("\nAfter cancelling 'r' and 'pi', the calculation simplifies to:")
    print(f"Percentage = ({num_arcs} / {angle_denominator}) / {divisor_coefficient} * {percentage_multiplier}")
    print(f"           = {num_arcs / angle_denominator} / {divisor_coefficient} * {percentage_multiplier}")
    print(f"           = {ratio} * {percentage_multiplier}")
    print(f"           = {final_percentage}")

# Execute the function to print the solution.
solve_cube_locus_problem()