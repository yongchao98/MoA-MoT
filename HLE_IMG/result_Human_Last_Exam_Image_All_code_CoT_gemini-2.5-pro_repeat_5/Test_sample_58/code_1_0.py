import math

def solve_cylinder_height():
    """
    This function explains the derivation for the height of the cylinder (h)
    in terms of its radius (r) and an angle (theta).
    """

    # The final formula is derived symbolically. We will print the steps.
    # The derivation assumes theta is in radians.

    print("Derivation for the height 'h' of the cylinder:")
    print("-" * 50)

    # Step 1: Calculate the length of the chord on the top face.
    # Let the vertex of the angle be P. Let the diameter be PA and the other chord be PB.
    # In the right triangle PAB (right-angled at B), the hypotenuse PA = 2*r.
    # The length of chord PB is PA * cos(theta).
    print("Step 1: Calculate the length of the chord on the top face (let's call it L_chord).")
    print("In the right triangle formed on the top circular face, the diameter (2*r) is the hypotenuse.")
    print("The chord's length is the adjacent side to the angle theta.")
    print("L_chord = (2 * r) * cos(theta)")
    print()

    # Step 2: Calculate the length of the geodesic path on the cylinder's surface.
    # This path connects the end of the chord (B) to the point on the bottom rim (Q).
    # Q is vertically below the end of the diameter (A).
    # When the cylinder is unrolled, this forms a right triangle.
    # The height of this triangle is h.
    # The base is the arc length between A and B.
    # The central angle for this arc is 2*theta.
    # The arc length is r * (2*theta).
    print("Step 2: Calculate the length of the path on the curved surface (let's call it L_geodesic).")
    print("When the cylinder surface is unrolled into a rectangle, this path becomes the hypotenuse of a right triangle.")
    print("The legs of this right triangle are the height 'h' and the arc length on the rim.")
    print("The arc length is calculated as r * (2 * theta).")
    print("Using the Pythagorean theorem:")
    print("L_geodesic^2 = h^2 + (2 * r * theta)^2")
    print()

    # Step 3: Equate the lengths and solve for h.
    # The key assumption is that L_chord = L_geodesic.
    print("Step 3: Assume the chord length equals the geodesic path length and solve for h.")
    print("L_chord^2 = L_geodesic^2")
    print("(2 * r * cos(theta))^2 = h^2 + (2 * r * theta)^2")
    print("4 * r^2 * cos(theta)^2 = h^2 + 4 * r^2 * theta^2")
    print("h^2 = 4 * r^2 * cos(theta)^2 - 4 * r^2 * theta^2")
    print("h^2 = 4 * r^2 * (cos(theta)^2 - theta^2)")
    print("h = sqrt(4 * r^2 * (cos(theta)^2 - theta^2))")
    print("-" * 50)

    # Step 4: Final Answer
    print("The final expression for the height h is:")
    # The numbers in the equation are 2 and 2 (as exponents).
    print("h = 2*r*sqrt(cos(theta)**2 - theta**2)")

solve_cylinder_height()