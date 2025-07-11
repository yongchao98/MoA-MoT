import numpy as np

def solve_fortress_problem_for_sphere():
    """
    Demonstrates that a finite number of guards on a sphere's surface
    cannot guard its entire exterior.

    This function sets up 4 guards in a regular tetrahedral formation on the
    surface of a unit sphere. It then calculates the position of a vertex
    of the polyhedron formed by the tangent planes at the guard locations.
    It shows this vertex is outside the sphere, and then finds a specific
    "blind spot" in the exterior that no guard can see.
    """
    print("Demonstration for k=4 guards in a regular tetrahedral arrangement.\n")

    # 1. Define the positions of 4 guards on the unit sphere.
    # These points form the vertices of a regular tetrahedron.
    g1 = np.array([np.sqrt(8/9), 0, -1/3])
    g2 = np.array([-np.sqrt(2/9), np.sqrt(6/9), -1/3])
    g3 = np.array([-np.sqrt(2/9), -np.sqrt(6/9), -1/3])
    g4 = np.array([0, 0, 1])
    
    guards = [g1, g2, g3, g4]

    print("Guard positions (vertices of an inscribed tetrahedron):")
    for i, g in enumerate(guards):
        print(f" g{i+1}: {np.round(g, 3)}, Magnitude: {np.linalg.norm(g):.1f}")
    
    print("-" * 30)

    # 2. Find a vertex of the polyhedron formed by the tangent planes.
    # A vertex is the intersection of 3 tangent planes. Let's find the
    # vertex v formed by the planes tangent at g1, g2, and g3.
    # The equations are <gi, v> = 1 for i=1,2,3.
    # This is a linear system A*v = b
    A = np.array([g1, g2, g3])
    b = np.array([1, 1, 1])
    
    vertex = np.linalg.solve(A, b)
    vertex_magnitude = np.linalg.norm(vertex)

    print(f"A vertex of the resulting 'unseen' polyhedron is at:")
    print(f" v = {np.round(vertex, 3)}")
    print(f"\nThe magnitude (distance from origin) of this vertex is:")
    print(f" |v| = {vertex_magnitude:.4f}")
    print("Since |v| > 1, this vertex is outside the unit ball.")

    print("-" * 30)

    # 3. Identify a blind spot.
    # Any point inside the polyhedron of unseen points is, by definition, unseen.
    # If we find such a point that is ALSO outside the ball, we have found a blind spot.
    # Let's test the point p, which is on the line from the origin to the vertex,
    # but still outside the ball. E.g., p = v * (2/3)
    
    blind_spot = vertex * (2 / 3.0)
    blind_spot_magnitude = np.linalg.norm(blind_spot)
    
    print(f"Checking for a blind spot at p = {np.round(blind_spot, 3)}")
    print(f"Magnitude of p is |p| = {blind_spot_magnitude:.4f}")
    print("Since |p| > 1, this point is in the area that needs to be guarded.\n")

    print("Checking if any guard can see point p:")
    can_see = False
    for i, g in enumerate(guards):
        # A guard g can see a point p if the dot product <g, p> >= 1.
        dot_product = np.dot(g, blind_spot)
        print(f" <g{i+1}, p> = {dot_product:.4f}")
        if dot_product >= 1:
            can_see = True

    print("\nConclusion:")
    if not can_see:
        print("No guard can see the point p.")
        print("This demonstrates that 4 guards are not sufficient.")
        print("This reasoning extends to any finite number of guards.")
    else:
        # This case should not be reached with the chosen point.
        print("An error occurred, the point was unexpectedly visible.")


solve_fortress_problem_for_sphere()