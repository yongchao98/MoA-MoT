def solve_geometry_ratio():
    """
    This script derives the ratio BM/MI in a triangle ABC.

    It explains the steps using geometric theorems and prints the final result.
    """

    # Define symbolic representations for the side lengths of the triangle ABC
    # a is the length of side BC (opposite to vertex A)
    # b is the length of side AC (opposite to vertex B)
    # c is the length of side AB (opposite to vertex C)
    side_a = 'a'
    side_b = 'b'
    side_c = 'c'

    print("Derivation of the ratio BM / MI:")
    print("-" * 40)

    # Step 1: Apply the Trillium Theorem
    print("Step 1: The Trillium Theorem")
    print("Let M be the point where the angle bisector of angle B intersects the circumcircle.")
    print("A key geometric property, the Trillium Theorem, states that M is equidistant from vertices A, C, and the incenter I.")
    print("This gives us the equality: MA = MC = MI.\n")

    # Step 2: Apply Ptolemy's Theorem
    print("Step 2: Ptolemy's Theorem")
    print("The points A, B, C, and M lie on the circumcircle, forming a cyclic quadrilateral ABCM.")
    print("Ptolemy's Theorem for this quadrilateral states that the product of the diagonals equals the sum of the products of the opposite sides:")
    print("AC * BM = AB * CM + BC * AM\n")

    # Step 3: Substitute known values
    print("Step 3: Substitute and Simplify")
    print(f"We substitute the side lengths (AC={side_b}, AB={side_c}, BC={side_a}) into the equation:")
    print(f"{side_b} * BM = {side_c} * CM + {side_a} * AM\n")

    print("Now, we use the result from the Trillium Theorem (MA = MC = MI) to replace AM and CM with MI:")
    print(f"{side_b} * BM = {side_c} * MI + {side_a} * MI\n")

    print("Factoring out MI from the right side gives:")
    print(f"{side_b} * BM = ({side_a} + {side_c}) * MI\n")

    # Step 4: Final Ratio
    print("Step 4: The Final Ratio")
    print("To find the desired ratio, we divide both sides by MI and by b:")
    print("(BM / MI) = (a + c) / b\n")

    # Final Answer Output
    print("-" * 40)
    print("The final expression for the ratio is:")
    # The prompt requires printing each component of the final equation
    numerator_part1 = side_a
    numerator_part2 = side_c
    denominator = side_b
    print(f"BM / MI = ({numerator_part1} + {numerator_part2}) / {denominator}")

# Execute the function to print the solution
solve_geometry_ratio()