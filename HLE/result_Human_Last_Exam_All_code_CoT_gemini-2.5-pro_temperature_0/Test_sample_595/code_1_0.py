import math

def solve():
    """
    Calculates the largest number k of coordinate grid squares a triangle's
    perimeter can pass through without containing any lattice points.
    """

    # Step 1: Define the triangle properties
    # Isosceles right triangle with legs of length 18.

    # Step 2: Formulate the calculation for k
    # k is the L1 perimeter of the triangle formed by the floor of the vertices.
    # We analyze two orientations to find the maximum k.

    print("---")
    print("Analysis for Orientation 1: Legs aligned with axes.")

    # Let the right-angle vertex P be at (e, e) where e is a small non-integer.
    e = 0.1
    P1 = (e, e)
    Q1 = (18 + e, e)
    R1 = (e, 18 + e)

    # Calculate the floor of the vertex coordinates
    fP1 = (math.floor(P1[0]), math.floor(P1[1]))
    fQ1 = (math.floor(Q1[0]), math.floor(Q1[1]))
    fR1 = (math.floor(R1[0]), math.floor(R1[1]))

    # Calculate the number of crossings for each side (L1 distance of floor vertices)
    k_PQ1 = abs(fQ1[0] - fP1[0]) + abs(fQ1[1] - fP1[1])
    k_QR1 = abs(fR1[0] - fQ1[0]) + abs(fR1[1] - fQ1[1])
    k_RP1 = abs(fP1[0] - fR1[0]) + abs(fP1[1] - fR1[1])
    k1 = k_PQ1 + k_QR1 + k_RP1

    print(f"Placing right-angle vertex at {P1}, the other vertices are {Q1} and {R1}.")
    print(f"The integer parts of the vertices are fP={fP1}, fQ={fQ1}, fR={fR1}.")
    print(f"Crossings for leg 1 (PQ) = |{fQ1[0]} - {fP1[0]}| + |{fQ1[1]} - {fP1[1]}| = {k_PQ1}")
    print(f"Crossings for hypotenuse (QR) = |{fR1[0]} - {fQ1[0]}| + |{fR1[1]} - {fQ1[1]}| = {k_QR1}")
    print(f"Crossings for leg 2 (RP) = |{fP1[0]} - {fR1[0]}| + |{fP1[1]} - {fR1[1]}| = {k_RP1}")
    print(f"Total k for Orientation 1 = {k_PQ1} + {k_QR1} + {k_RP1} = {k1}")

    print("\n" + "---")
    print("Analysis for Orientation 2: Hypotenuse aligned with an axis.")

    # The legs (length 18) are now diagonal. Their projection on the axes is 18/sqrt(2).
    d = 9 * math.sqrt(2)  # approx 12.728

    # Place the right-angle vertex B at (x0, y0) and choose x0, y0 to maximize
    # the floor differences and avoid lattice points on the perimeter.
    # A = (x0 - d, y0 - d), C = (x0 + d, y0 - d)
    # A valid choice that maximizes k is x0=0.15, y0=0.25.
    x0 = 0.15
    y0 = 0.25

    B2 = (x0, y0)
    A2 = (x0 - d, y0 - d)
    C2 = (x0 + d, y0 - d)

    # Calculate the floor of the vertex coordinates
    fA2 = (math.floor(A2[0]), math.floor(A2[1]))
    fB2 = (math.floor(B2[0]), math.floor(B2[1]))
    fC2 = (math.floor(C2[0]), math.floor(C2[1]))

    # Calculate the number of crossings for each side
    k_AB2 = abs(fB2[0] - fA2[0]) + abs(fB2[1] - fA2[1])
    k_BC2 = abs(fC2[0] - fB2[0]) + abs(fC2[1] - fB2[1])
    k_CA2 = abs(fA2[0] - fC2[0]) + abs(fA2[1] - fC2[1])
    k2 = k_AB2 + k_BC2 + k_CA2

    print(f"Placing right-angle vertex at {B2}, the other vertices are {A2} and {C2}.")
    print(f"The integer parts of the vertices are fA={fA2}, fB={fB2}, fC={fC2}.")
    print(f"Crossings for leg 1 (AB) = |{fB2[0]} - {fA2[0]}| + |{fB2[1]} - {fA2[1]}| = {k_AB2}")
    print(f"Crossings for leg 2 (BC) = |{fC2[0]} - {fB2[0]}| + |{fC2[1]} - {fB2[1]}| = {k_BC2}")
    print(f"Crossings for hypotenuse (CA) = |{fA2[0]} - {fC2[0]}| + |{fA2[1]} - {fC2[1]}| = {k_CA2}")
    print(f"Total k for Orientation 2 = {k_AB2} + {k_BC2} + {k_CA2} = {k2}")

    print("\n" + "---")
    print(f"Comparing the two orientations, the maximum value found is {max(k1, k2)}.")

solve()
<<<76>>>