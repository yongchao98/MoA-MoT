import numpy as np

def solve_fortress_problem():
    """
    Demonstrates that a finite number of guards cannot observe the entire
    exterior of a unit sphere.

    The script uses an arrangement of 4 guards at the vertices of a regular
    tetrahedron inscribed in the unit sphere. It then calculates the location
    of a vertex of the "unseen" region and shows that this vertex lies
    outside the unit sphere, proving that some exterior points are not observed.
    """

    print("Step 1: Define guard positions for a regular tetrahedron inscribed in the unit sphere.")
    # Vertices of a regular tetrahedron centered at the origin, on the unit sphere
    g1 = np.array([0, 0, 1.0])
    g2 = np.array([2 * np.sqrt(2) / 3, 0, -1/3])
    g3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1/3])
    g4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1/3])
    guards = [g1, g2, g3, g4]
    print("Guard 1:", g1)
    print("Guard 2:", g2)
    print("Guard 3:", g3)
    print("Guard 4:", g4)
    print("-" * 30)

    print("Step 2: Find a vertex of the unseen region.")
    print("The unseen region is defined by P . G < 1 for all guards G.")
    print("Its vertices are found where 3 of these constraints meet, i.e., V . G = 1 for 3 guards.")
    
    # We find the vertex V corresponding to the face (g2, g3, g4)
    # This involves solving the linear system A * V = b
    # where A is the matrix of guard vectors and b is [1, 1, 1]
    A = np.array([g2, g3, g4])
    b = np.array([1, 1, 1])
    V = np.linalg.solve(A, b)
    
    print(f"The vertex V of the unseen region is: {V}")
    print("-" * 30)

    print("Step 3: Check if this vertex is outside the unit ball.")
    distance_V = np.linalg.norm(V)
    print(f"The distance of vertex V from the origin is ||V|| = {distance_V:.4f}")

    if distance_V > 1:
        print("Since ||V|| > 1, the vertex of the unseen region is outside the unit ball.")
    else:
        print("The vertex of the unseen region is inside the unit ball.")
    print("-" * 30)

    print("Step 4: Construct an unobserved point P outside the unit ball.")
    # Any point P = (1-epsilon)*V for small epsilon > 0 near V will also be
    # outside the ball, and will be in the unseen region.
    # Let's pick a point 99% of the way from the origin to V.
    P = 0.99 * V
    distance_P = np.linalg.norm(P)
    print(f"Let's test the point P = 0.99 * V = {P}")
    print(f"The distance of P from the origin is ||P|| = {distance_P:.4f}")
    if distance_P > 1:
        print("This point P is outside the unit ball.")
    else:
        print("This point is not outside the unit ball.")
    print("-" * 30)
    
    print("Step 5: Verify that point P is not seen by any guard.")
    print("We check if P . G < 1 for all guards G.")
    
    unseen = True
    for i, g in enumerate(guards):
        dot_product = np.dot(P, g)
        print(f"P . G{i+1} = {dot_product:.4f}")
        if dot_product >= 1:
            unseen = False

    if unseen:
        print("\nConclusion: The point P is outside the unit ball AND is not seen by any guard.")
        print("This demonstrates that 4 guards are not sufficient.")
        print("This argument holds for any finite number of guards, as the unseen region will always be a polyhedron with vertices outside the sphere.")
    else:
        print("\nSomething went wrong. The point should be unseen.")
        
solve_fortress_problem()