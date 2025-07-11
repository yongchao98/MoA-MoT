def solve_geometry_ratio():
    """
    This function solves the geometry problem and prints the final expression for the ratio.

    The derivation involves the following key steps:
    1.  Using the Trillium Theorem, which states that the intersection M of the angle
        bisector BI with the circumcircle is the center of the circle through A, I, and C.
        This gives the crucial equality MA = MI = MC.
    
    2.  The ratio to be found, BM/MI, can therefore be rewritten as BM/MC.
    
    3.  Applying the Law of Sines to the triangle BMC, we get:
        BM / sin(∠BCM) = MC / sin(∠CBM)
        This implies BM/MC = sin(∠BCM) / sin(∠CBM).

    4.  The angles are identified using properties of cyclic quadrilaterals and angle bisectors:
        - ∠CBM = B/2
        - ∠BCM = C + B/2

    5.  The ratio simplifies to sin(C + B/2) / sin(B/2).

    6.  This trigonometric expression is proven to be equivalent to (a+c)/b by using
        the sum-to-product and half-angle identities along with the Law of Sines on triangle ABC.
    """

    # Representing the side lengths symbolically
    side_a = "a"
    side_b = "b"
    side_c = "c"
    
    # Printing the final expression for the ratio BM/MI
    print("The ratio BM/MI in terms of the side lengths a, b, and c is:")
    print(f"({side_a} + {side_c}) / {side_b}")

if __name__ == "__main__":
    solve_geometry_ratio()