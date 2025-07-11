import sys

def solve_geometry_ratio():
    """
    This function provides a step-by-step derivation of the ratio BM/MI
    in triangle ABC and prints the final result.
    """

    print("Problem: In a triangle ABC, the angle bisectors of angle BAC and angle ABC intersect at I.")
    print("Let M be the intersection of the line BI with the circumcircle of triangle ABC.")
    print("Express the ratio BM/MI in terms of the side lengths a, b, and c.\n")

    print("--- Step 1: The Key Property of Point M ---")
    print("A crucial geometric property relates the point M to the incenter I and vertices A and C.")
    print("This property is that M is equidistant from A, I, and C. That is: MA = MI = MC.")
    print("This is sometimes known as the 'Fact 5' lemma or the 'center of the circle through A, I, C'.")
    print("Let's prove this briefly:")
    print("  - Let the angles at vertices A, B, C be alpha, beta, gamma.")
    print("  - Since BI is an angle bisector, angle CBM = beta/2.")
    print("  - In the circumcircle, arc MC subtends both angle MAC and angle CBM. Thus, angle MAC = beta/2.")
    print("  - The angle AIM is the exterior angle of triangle ABI, so AIM = angle IAB + angle IBA = alpha/2 + beta/2.")
    print("  - The angle IAM = angle IAC + angle CAM = alpha/2 + beta/2.")
    print("  - Since angle AIM = angle IAM, triangle AIM is isosceles with MA = MI.")
    print("  - A similar proof shows that triangle CIM is also isosceles with MC = MI.")
    print("  - Therefore, we have proven that MA = MI = MC.\n")

    print("--- Step 2: Applying Ptolemy's Theorem ---")
    print("The points A, B, C, and M all lie on the circumcircle, forming a cyclic quadrilateral ABCM.")
    print("Ptolemy's Theorem for a cyclic quadrilateral states:")
    print("(product of diagonals) = (sum of products of opposite sides)")
    print("For quadrilateral ABCM, the diagonals are AC and BM. The pairs of opposite sides are (AB, CM) and (BC, AM).")
    print("The equation is: AC * BM = AB * CM + BC * AM\n")

    print("--- Step 3: Substituting and Solving ---")
    print("Let the side lengths opposite to A, B, C be a, b, c respectively.")
    print("So, AB = c, BC = a, and AC = b.")
    print("Substituting these into Ptolemy's equation:")
    print("b * BM = c * CM + a * AM")
    print("Now, we use our key property from Step 1: MA = MC = MI.")
    print("Substitute MI for both AM and CM in the equation:")
    print("b * BM = c * MI + a * MI")
    print("Factor out MI on the right side:")
    print("b * BM = (a + c) * MI")
    print("Finally, rearrange to find the desired ratio BM / MI:")
    print("BM / MI = (a + c) / b\n")

    print("--- Final Expression ---")
    side_a = 'a'
    side_b = 'b'
    side_c = 'c'
    # The final equation consists of these parts.
    # We print them to form the final expression as requested.
    print(f"The ratio BM/MI is equal to ({side_a} + {side_c}) / {side_b}")

if __name__ == '__main__':
    solve_geometry_ratio()