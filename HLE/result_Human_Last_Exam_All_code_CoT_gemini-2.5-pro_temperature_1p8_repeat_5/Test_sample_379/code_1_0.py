import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus on a cube's surface and a related percentage.
    """
    print("This problem asks for the length of a specific locus of points on the surface of a cube.")
    print("My plan is to calculate this length by 'unfolding' the cube's surface for each relevant face.")
    print("-" * 20)

    # Let r be the side length. The calculations are proportional to r and pi.
    # We can represent lengths as a coefficient multiplied by (pi * r).

    # Length on the Front face (adjacent to P's edge)
    # The unfolding of the Front and Right faces shows an arc on the Front face.
    # Calculation shows this arc subtends an angle of pi/3.
    # Arc Length = radius * angle = r * (pi/3)
    len_coeff_F = 1/3

    # Length on the Right face (adjacent to P's edge)
    # This is symmetrical to the Front face's arc.
    len_coeff_R = 1/3

    # Length on the Top face
    # The locus on the Top face consists of two arcs.
    # One via the Front-Top unfolding (length r*pi/3).
    # One via the Right-Top unfolding (length r*pi/3).
    # Total length on Top = r*pi/3 + r*pi/3 = 2*pi*r/3.
    len_coeff_T = 2/3

    # Length on the Bottom face (by symmetry with Top)
    len_coeff_D = 2/3

    # The far faces (Left and Back) are at a minimum distance of r from P.
    # The locus only touches them at single points, which have zero length.

    # --- Total Length Calculation ---
    print("The total length of the locus C is the sum of lengths on each contributing face.")
    print("Let L_F, L_R, L_T, L_D be the lengths on the Front, Right, Top, and Bottom faces, and r be the side length.")
    print(f"L_F = ({len_coeff_F:.2f}...) * pi * r = (1/3) * pi * r")
    print(f"L_R = ({len_coeff_R:.2f}...) * pi * r = (1/3) * pi * r")
    print(f"L_T = ({len_coeff_T:.2f}...) * pi * r = (2/3) * pi * r")
    print(f"L_D = ({len_coeff_D:.2f}...) * pi * r = (2/3) * pi * r")

    # Summing the coefficients for the final equation
    total_coeff = len_coeff_F + len_coeff_R + len_coeff_T + len_coeff_D
    
    print("\nThe total length C is the sum:")
    # We display each number in the final equation as requested.
    print(f"C = (1/3)*pi*r + (1/3)*pi*r + (2/3)*pi*r + (2/3)*pi*r")
    print(f"C = ({len_coeff_F:.2f} + {len_coeff_R:.2f} + {len_coeff_T:.2f} + {len_coeff_D:.2f}) * pi * r")
    print(f"C = ({total_coeff:.2f}) * pi * r")
    print("Since (1/3) + (1/3) + (2/3) + (2/3) = 6/3 = 2, we have:")
    print("C = 2 * pi * r")

    # --- Final Ratio Calculation ---
    # The problem asks to divide C by 2*pi*r and give the percentage.
    # Let r=1 for the calculation, as it cancels out.
    r = 1.0
    total_length_C = total_coeff * math.pi * r
    divisor = 2 * math.pi * r
    ratio = total_length_C / divisor
    percentage = int(round(ratio * 100))

    print("\nNow, we divide the length of C by 2*pi*r:")
    print(f"Ratio = (2 * pi * r) / (2 * pi * r) = {ratio:.0f}")

    print("\nFinally, we express this ratio as a whole number percentage:")
    print(f"Percentage = {ratio:.0f} * 100% = {percentage}%")

solve_cube_locus_problem()