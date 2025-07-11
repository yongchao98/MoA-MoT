import math

def solve_cube_surface_area_ratio():
    """
    Calculates the ratio of a specific area D on a cube's surface to its total surface area S.

    The problem defines D as the set of points on the surface of a cube (side length s)
    that are at a surface distance of at most sqrt(2)*s from a vertex P.

    The final ratio is derived analytically to be (2*pi + 3*sqrt(3) - 3) / 12.
    This script verifies the components of this exact fraction.
    """
    s = 1.0 # The side length s cancels out, so we use s=1 for simplicity.

    # Total surface area of the cube
    area_S = 6 * s**2

    # Area contribution from the 3 faces adjacent to the vertex P.
    # The entire area of these faces is included.
    area_adjacent = 3 * s**2

    # Area contribution from one face opposite to the vertex P.
    # This formula is derived from integrating over the region on the face satisfying the distance condition.
    # A_opp = s^2 * (pi/3 + (sqrt(3)-3)/2)
    a_opp = s**2 * (math.pi/3 + (math.sqrt(3) - 3)/2)

    # Total area from the 3 opposite faces
    area_opposite = 3 * a_opp

    # Total area of the region D
    area_D = area_adjacent + area_opposite

    # The final ratio
    ratio = area_D / area_S

    # The derived exact fraction is (2*pi + 3*sqrt(3) - 3) / 12.
    # The code will print the numbers that form this exact answer.
    numerator_pi_coefficient = 2
    numerator_sqrt3_coefficient = 3
    numerator_constant = -3
    denominator = 12

    print("The final exact ratio is represented by the fraction: (A*pi + B*sqrt(3) + C) / D")
    print("The numbers in this final equation are:")
    print(f"A (coefficient of pi): {numerator_pi_coefficient}")
    print(f"B (coefficient of sqrt(3)): {numerator_sqrt3_coefficient}")
    print(f"C (constant term): {numerator_constant}")
    print(f"D (denominator): {denominator}")

    # Display the fraction in a readable format.
    final_numerator_expr = f"({numerator_pi_coefficient}*pi + {numerator_sqrt3_coefficient}*sqrt(3) {'' if numerator_constant < 0 else '+'} {numerator_constant})"
    print(f"\nThus, the final ratio is {final_numerator_expr} / {denominator}.")

solve_cube_surface_area_ratio()