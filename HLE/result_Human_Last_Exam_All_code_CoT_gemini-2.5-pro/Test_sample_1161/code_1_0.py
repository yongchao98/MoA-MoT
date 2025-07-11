import numpy as np

def solve_fortress_problem_sphere():
    """
    Demonstrates that a finite number of guards is insufficient to observe
    the exterior of a sphere.
    """
    # 1. Define the guard positions.
    # We place 4 guards at the vertices of a regular tetrahedron
    # inscribed in the unit sphere.
    g1 = np.array([np.sqrt(8)/3, 0, -1/3])
    g2 = np.array([-np.sqrt(2)/3, np.sqrt(6)/3, -1/3])
    g3 = np.array([-np.sqrt(2)/3, -np.sqrt(6)/3, -1/3])
    g4 = np.array([0, 0, 1])
    guards = [g1, g2, g3, g4]

    # 2. Define a test point `p` to check for visibility.
    # The tangent planes at g1, g2, g3 intersect at the point (0, 0, -3).
    # This point is a vertex of the tetrahedron that circumscribes the sphere.
    # Let's pick a point on the line from the origin to this vertex, but
    # still outside the unit ball.
    # The point p = (0, 0, -1.5) is outside the unit ball because its
    # distance from the origin is 1.5.
    p_unseen = np.array([0, 0, -1.5])
    
    print(f"Testing if the point p = {p_unseen.tolist()} is visible.")
    print(f"The magnitude of p is {np.linalg.norm(p_unseen)}, which is > 1, so it is outside the unit ball.\n")
    print("A guard at 'g' sees point 'p' if the dot product (p . g) >= 1.")
    print("-" * 60)

    is_seen = False
    for i, g in enumerate(guards):
        dot_product = np.dot(p_unseen, g)
        
        # The final "equation" is the check for visibility
        is_seen_by_guard = dot_product >= 1
        
        print(f"Guard {i+1} at {np.round(g, 2).tolist()}:")
        print(f"  Calculation: p . g = {dot_product:.4f}")
        # Output the numbers in the final check
        print(f"  Check: {dot_product:.4f} >= 1 is {is_seen_by_guard}")
        if is_seen_by_guard:
            is_seen = True
        print("-" * 60)

    print("\n--- Conclusion ---")
    if not is_seen:
        print("The point p is NOT seen by any of the 4 guards.")
        print("This demonstrates that 4 guards are not sufficient.")
    else:
        print("The point p was seen by at least one guard.")
        print("(This would indicate an error in the choice of p).")
        
    print("\nThis principle holds for any finite number of guards. A finite set of guards")
    print("defines a circumscribing polyhedron, which will always contain points")
    print("outside the sphere that remain unobserved.")

solve_fortress_problem_sphere()