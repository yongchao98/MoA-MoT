def solve_geometry_ratio():
    """
    This function prints the step-by-step derivation for the geometric ratio BM/MI.
    """

    # Symbolic representation of side lengths
    a = 'a'
    b = 'b'
    c = 'c'

    print("Problem: In a triangle ABC, the angle bisectors of angle BAC and angle ABC intersect at a point I.")
    print("Let M be the point where the angle bisector BI intersects the circumcircle of triangle ABC.")
    print("Express the ratio BM / MI in terms of the side lengths a, b, and c.")
    print("\n" + "="*50 + "\n")

    print("Step 1: Establishing the key equality MI = AM = CM (Trillium Theorem).")
    print("-" * 20)
    print("First, we prove that AM = CM.")
    print("Since BI is the angle bisector of angle B, it divides the arc AC of the circumcircle into two equal arcs.")
    print("Chords subtending equal arcs are equal, therefore AM = CM.")
    print("\nNext, we prove that AM = MI by showing triangle AIM is isosceles.")
    print("Let the angles of triangle ABC be A, B, and C.")
    print("  - Angle IAM = Angle IAC + Angle CAM = A/2 + Angle CBM = A/2 + B/2.")
    print("  - Angle AIM is the exterior angle of triangle AIB, so Angle AIM = Angle IAB + Angle IBA = A/2 + B/2.")
    print("Since Angle IAM = Angle AIM, triangle AIM is isosceles with AM = MI.")
    print("Thus, we have the crucial relationship: MI = AM = CM.")
    print("\n" + "="*50 + "\n")

    print("Step 2: Applying Ptolemy's Theorem.")
    print("-" * 20)
    print("The points A, B, C, and M form a cyclic quadrilateral ABCM.")
    print("Ptolemy's Theorem states: (AB * CM) + (BC * AM) = (AC * BM).")
    print(f"Substituting the side lengths (AB=c, BC=a, AC=b), we get:")
    print(f"({c} * CM) + ({a} * AM) = ({b} * BM)")
    print("\n" + "="*50 + "\n")

    print("Step 3: Solving for the ratio BM / MI.")
    print("-" * 20)
    print("From Step 1, we know AM = CM. Substitute CM with AM in the equation:")
    print(f"({c} * AM) + ({a} * AM) = ({b} * BM)")
    print("Factoring out AM:")
    print(f"({a} + {c}) * AM = {b} * BM")
    print("\nRearranging to find the ratio BM / AM:")
    print(f"BM / AM = ({a} + {c}) / {b}")
    print("\nSince MI = AM (from Step 1), we can substitute AM with MI to get the final answer.")
    print("The final expression for the ratio is:")
    print(f"BM / MI = ({a} + {c}) / {b}")
    print("\n" + "="*50 + "\n")

solve_geometry_ratio()