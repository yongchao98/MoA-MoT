import math

def calculate_locus_percentage():
    """
    This script solves the cube locus problem.

    The problem is to find the length of a curve C on a cube's surface and
    express a specific ratio as a whole number percentage.

    Through geometric analysis by unfolding the cube's faces, it can be shown
    that the total length of the curve C is exactly 2 * pi * r, where r is
    the side length of the cube.
    """

    # The coefficient for the total length of C is 2, for the expression 2 * pi * r.
    numerator_coefficient = 2

    # The problem asks to divide the length of C by the expression 2 * pi * r.
    # The coefficient for the divisor is also 2.
    denominator_coefficient = 2

    # The variables 'pi' and 'r' cancel out in the division.
    ratio = numerator_coefficient / denominator_coefficient

    # Convert the ratio to a whole number percentage.
    percentage = int(ratio * 100)

    # Print the final equation with all its numeric components, as requested.
    # The 'r' is included to show the full context of the original expressions.
    print(f"The total length of the locus C is {numerator_coefficient} * pi * r.")
    print(f"The problem asks to calculate: (Length of C) / ({denominator_coefficient} * pi * r) * 100%")
    print("\nDisplaying the final equation with the calculated numbers:")
    print(f"({numerator_coefficient} * pi * r) / ({denominator_coefficient} * pi * r) * 100% = {percentage}%")


calculate_locus_percentage()
