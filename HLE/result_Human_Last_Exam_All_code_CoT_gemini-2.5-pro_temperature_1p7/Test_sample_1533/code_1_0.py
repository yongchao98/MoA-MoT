import sympy

def solve_geometry_ratio():
    """
    This function provides a step-by-step derivation for the ratio BM/MI
    and prints the final result in terms of side lengths a, b, and c.
    """

    # Define symbols for clarity in explanation
    a, b, c = sympy.symbols('a b c')
    A, B, C = sympy.symbols('A B C') # Angles of triangle ABC

    print("Step 1: Relate the length of MI to another segment.")
    print("Let I be the incenter and M be the point where the angle bisector BI intersects the circumcircle.")
    print("A key property is that MI = MA. This is because triangle AIM is isosceles.")
    print("Proof for MI = MA:")
    print("  - Angle IAM = Angle IAC + Angle CAM = A/2 + Angle CBM = A/2 + B/2.")
    print("  - Angle AIM is the exterior angle of triangle AIB, so Angle AIM = Angle IAB + Angle IBA = A/2 + B/2.")
    print("  - Since Angle IAM = Angle AIM, triangle AIM is isosceles with MI = MA.")
    print("-" * 20)

    print("Step 2: Rewrite the ratio using the property from Step 1.")
    print("The required ratio BM / MI can be rewritten as BM / MA.")
    print("-" * 20)

    print("Step 3: Apply the Law of Sines to triangle ABM.")
    print("The angles of triangle ABM are:")
    print("  - Angle ABM = B/2 (since BM is the angle bisector).")
    print("  - Angle AMB = Angle ACB = C (angles subtended by the same arc AB).")
    print("  - Angle MAB = 180 - (B/2 + C) = A + B/2.")
    print("By the Law of Sines in triangle ABM:")
    print("  BM / sin(MAB) = MA / sin(ABM)")
    print("  BM / MA = sin(A + B/2) / sin(B/2)")
    print("-" * 20)

    print("Step 4: Convert the trigonometric expression into side lengths.")
    print("We can prove that sin(A + B/2) / sin(B/2) = (a + c) / b.")
    print("  (a + c) / b = (sin(A) + sin(C)) / sin(B)")
    print("  Using sum-to-product and double-angle identities, this simplifies to cos((A-C)/2) / sin(B/2).")
    print("  Similarly, sin(A + B/2) can be shown to be equal to cos((A-C)/2).")
    print("  Therefore, sin(A + B/2) / sin(B/2) = cos((A-C)/2) / sin(B/2).")
    print("This confirms that BM / MA = (a + c) / b.")
    print("-" * 20)

    print("Final Result:")
    print("The ratio BM / MI in terms of side lengths a, b, and c is:")
    
    # In the expression (a+c)/b, the coefficient for each variable is 1.
    print(f"BM / MI = ({1} * a + {1} * c) / ({1} * b)")
    print("\nOr more simply:")
    print("BM / MI = (a + c) / b")


solve_geometry_ratio()