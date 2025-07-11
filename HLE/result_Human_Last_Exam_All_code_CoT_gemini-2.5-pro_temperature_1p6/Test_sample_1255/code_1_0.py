import sympy

def check_p1_multiplication():
    """
    Demonstrates that the multiplication map for P^1 is not a morphism.

    In projective geometry, P^1 can be represented by homogeneous coordinates [x_0 : x_1],
    where not both are zero. The map from an affine coordinate x to homogeneous
    coordinates is x -> [1 : x]. The point at infinity is [0 : 1].

    The multiplication map m(x, y) = xy on the affine parts corresponds to a map on
    homogeneous coordinates. Let P = [x_0 : x_1] and Q = [y_0 : y_1].
    The product is m(P, Q) = [x_0*y_0 : x_1*y_1].

    A map is well-defined if the resulting coordinates are not both zero.
    We check for points of indeterminacy.
    """

    x0, x1, y0, y1 = sympy.symbols('x0 x1 y0 y1')

    # Coordinates of the product
    z0 = x0 * y0
    z1 = x1 * y1

    print("Investigating the multiplication map m([x0:x1], [y0:y1]) = [x0*y0 : x1*y1] on P^1 x P^1.")
    print(f"The resulting coordinates are z0 = {z0} and z1 = {z1}.")
    print("The map is undefined if both z0 and z1 are simultaneously zero.\n")

    # We need to solve the system of equations z0 = 0 and z1 = 0,
    # keeping in mind that for a point on P^1, (x0, x1) != (0, 0) and (y0, y1) != (0, 0).

    # Case 1: Point P1 is (0, inf), corresponding to affine x=0, affine y=inf
    # Affine x=0 corresponds to homogeneous coords [1 : 0]. So x0=1, x1=0.
    # Affine y=inf corresponds to homogeneous coords [0 : 1]. So y0=0, y1=1.
    P1_x0, P1_x1 = 1, 0
    P1_y0, P1_y1 = 0, 1
    res_z0 = z0.subs({x0: P1_x0, x1: P1_x1, y0: P1_y0, y1: P1_y1})
    res_z1 = z1.subs({x0: P1_x0, x1: P1_x1, y0: P1_y0, y1: P1_y1})
    print(f"Checking the point (x, y) = (0, infinity):")
    print(f"Homogeneous coords: P = [{P1_x0}:{P1_x1}], Q = [{P1_y0}:{P1_y1}]")
    print(f"Resulting coords: z0 = {res_z0}, z1 = {res_z1}")
    if res_z0 == 0 and res_z1 == 0:
        print("Result is [0:0], which is undefined. This is a point of indeterminacy.\n")

    # Case 2: Point P2 is (inf, 0), corresponding to affine x=inf, affine y=0
    # Affine x=inf corresponds to homogeneous coords [0 : 1]. So x0=0, x1=1.
    # Affine y=0 corresponds to homogeneous coords [1 : 0]. So y0=1, y1=0.
    P2_x0, P2_x1 = 0, 1
    P2_y0, P2_y1 = 1, 0
    res_z0 = z0.subs({x0: P2_x0, x1: P2_x1, y0: P2_y0, y1: P2_y1})
    res_z1 = z1.subs({x0: P2_x0, x1: P2_x1, y0: P2_y0, y1: P2_y1})
    print(f"Checking the point (x, y) = (infinity, 0):")
    print(f"Homogeneous coords: P = [{P2_x0}:{P2_x1}], Q = [{P2_y0}:{P2_y1}]")
    print(f"Resulting coords: z0 = {res_z0}, z1 = {res_z1}")
    if res_z0 == 0 and res_z1 == 0:
        print("Result is [0:0], which is undefined. This is a point of indeterminacy.\n")
    
    print("Conclusion: Because the multiplication map is not a well-defined morphism everywhere on P^1 x P^1,")
    print("the underlying scheme P^1 of the logarithmic multiplicative group is not a group scheme.")
    print("This serves as a counterexample to the statement.")


if __name__ == '__main__':
    check_p1_multiplication()