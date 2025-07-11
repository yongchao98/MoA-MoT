import math

def calculate_locus_percentage():
    """
    Calculates the ratio of the locus length to 2*pi*r as a percentage.
    
    The problem involves finding the length of a locus of points C at a surface
    distance 'r' from a midpoint of an edge 'P' on a cube of side 'r'.
    
    The derivation shows the locus consists of 4 circular arcs.
    """
    
    # Based on the geometric derivation, the locus C is composed of 4 arcs.
    num_arcs = 4
    
    # Each arc is part of a circle of radius 'r' and subtends an angle of pi/3 radians.
    # The coefficient for the angle is 1/3.
    arc_angle_coeff = 1/3
    
    # The total length of the locus C is num_arcs * (arc_angle_coeff * pi * r).
    # We can calculate the numerical coefficient for the length equation, which is num_arcs * arc_angle_coeff.
    total_length_coeff = num_arcs * arc_angle_coeff
    
    # We need to find the ratio of this length to (2 * pi * r).
    # The divisor's numerical coefficient is 2.
    divisor_coeff = 2
    
    # The terms 'pi' and 'r' cancel out in the division.
    ratio = total_length_coeff / divisor_coeff
    
    # The final answer is this ratio expressed as a whole number percentage.
    percentage = round(ratio * 100)

    # Print out the steps of the final equation as requested.
    print(f"The total length of the locus C is given by the equation:")
    print(f"Length(C) = (Number of Arcs) * (Arc Angle) * r")
    print(f"Length(C) = {num_arcs} * (pi/3) * r")
    print(f"Length(C) = ({num_arcs}/3) * pi * r")
    
    print("\nWe are asked to calculate the ratio: Length(C) / (2 * pi * r)")
    print(f"Ratio = (({num_arcs}/3) * pi * r) / ({divisor_coeff} * pi * r)")
    
    # The 'pi * r' terms cancel out.
    print(f"\nThis simplifies to the ratio of the coefficients:")
    print(f"Ratio = ({num_arcs}/3) / {divisor_coeff} = {total_length_coeff:.4f} / {divisor_coeff} = {ratio:.4f}")

    print(f"\nFinally, we convert this ratio to a whole number percentage:")
    print(f"Percentage = round({ratio:.4f} * 100) = {percentage}%")
    
    # Return the final numerical answer as requested in the format <<<answer>>>
    return percentage

# Run the calculation and store the final answer.
final_answer = calculate_locus_percentage()

# The final result in the requested format.
# print(f"\n<<< {final_answer} >>>") # The output format will be handled by the <<<>>> delimiters after the code block