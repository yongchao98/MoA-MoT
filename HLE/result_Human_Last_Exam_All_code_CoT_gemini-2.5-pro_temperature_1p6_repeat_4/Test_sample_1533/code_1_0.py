def solve_geometry_ratio():
    """
    This function explains the derivation and prints the final expression for the ratio BM/MI.
    """

    # Symbolic representations of the side lengths of triangle ABC
    # a is the length of side BC (opposite to vertex A)
    # b is the length of side AC (opposite to vertex B)
    # c is the length of side AB (opposite to vertex C)
    side_a = 'a'
    side_b = 'b'
    side_c = 'c'

    # The problem is to express the ratio BM/MI in terms of a, b, and c.

    # Step 1: Prove that MA = MI = MC.
    # Let A, B, C denote the angles of the triangle at vertices A, B, C.
    # I is the incenter, so AI and BI are angle bisectors.
    # M is on the circumcircle and on the line BI.
    #
    # Proof that MA = MI:
    # In triangle AIM, we analyze its angles.
    # angle MAI = angle MAC + angle CAI
    # M is on the circumcircle, so angle MAC and angle MBC subtend the same arc MC.
    # Thus, angle MAC = angle MBC.
    # BI is the angle bisector of angle B, so angle MBC = B/2.
    # AI is the angle bisector of angle A, so angle CAI = A/2.
    # Therefore, angle MAI = B/2 + A/2.
    #
    # Now consider angle MIA. It's an exterior angle to triangle AIB.
    # angle MIA = angle IAB + angle IBA
    # angle IAB = A/2 and angle IBA = B/2.
    # Therefore, angle MIA = A/2 + B/2.
    # Since angle MAI = angle MIA, triangle AIM is an isosceles triangle with MA = MI.
    #
    # Proof that MA = MC:
    # Since the line BM is the angle bisector of angle ABC, it bisects the arc AC of the circumcircle
    # that does not contain the vertex B. The chords subtending these equal arcs must be equal.
    # Thus, MA = MC.
    # Combining the results, we have established that MA = MC = MI.

    # Step 2: Apply Ptolemy's Theorem.
    # The vertices A, B, C, and M lie on the circumcircle, forming a cyclic quadrilateral ABCM.
    # Ptolemy's theorem states: AB * CM + BC * AM = AC * BM.

    # Step 3: Substitute and solve.
    # Using our side length notation: c * CM + a * AM = b * BM.
    # From Step 1, we substitute MI for both AM and CM:
    # c * MI + a * MI = b * BM
    # (a + c) * MI = b * BM

    # Step 4: Isolate the desired ratio.
    # Dividing both sides by MI and by b, we get:
    # BM / MI = (a + c) / b

    # Final result presentation
    numerator = f"{side_a} + {side_c}"
    denominator = f"{side_b}"

    print("The ratio of the lengths of the segments BM and MI is derived using geometric properties.")
    print("The final expression in terms of the side lengths a, b, and c is:")
    print(f"BM / MI = ({numerator}) / {denominator}")

# Execute the function to get the solution
solve_geometry_ratio()