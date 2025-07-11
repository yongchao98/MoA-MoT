def solve_icosahedron_problem():
    """
    This script explains the reasoning to determine the shape of the water
    surface in a half-filled icosahedron tank standing on a face.
    """

    print("--- Problem Analysis ---")
    print("We need to find the shape of the water surface in an icosahedron tank.")
    print("Key conditions: The tank is an icosahedron, it stands on one face (base is a horizontal triangle), and it's half-filled.")
    print("\n")

    print("--- Step 1: The 'Half-Filled' Condition ---")
    print("An icosahedron is a centrally symmetric solid.")
    print("This means it has a geometric center, and for every point on the shape, there is a corresponding point directly opposite through the center.")
    print("When such a shape is half-filled, the surface of the water forms a horizontal plane passing exactly through this geometric center.")
    print("\n")

    print("--- Step 2: Analyzing the Symmetries of the Cross-Section ---")
    print("The shape of the water surface is the shape of the cross-section created by this central plane.")
    print("Let's consider the symmetries of the icosahedron when it stands on a face:")
    print("  1. Rotational Symmetry: The vertical axis through the center of the triangular base is an axis of 3-fold rotational symmetry. Since the water surface is a horizontal plane perpendicular to this axis, the shape of the surface must also have 3-fold rotational symmetry.")
    print("  2. Point Symmetry: The water surface plane passes through the icosahedron's center of symmetry. Therefore, the cross-sectional shape must also be centrally symmetric (have point symmetry).")
    print("\n")

    print("--- Step 3: Deducing the Final Shape ---")
    print("The plane cuts 6 of the icosahedron's edges around its 'equator'. The intersection points form the vertices of the water's surface shape.")
    print("A shape with 6 vertices that has both 3-fold rotational symmetry and point symmetry must be a regular hexagon.")
    print("Therefore, all its sides are of equal length and all its internal angles are equal.")
    print("\n")

    print("--- Conclusion ---")
    print("The shape of the water surface is a REGULAR HEXAGON.")


if __name__ == "__main__":
    solve_icosahedron_problem()
