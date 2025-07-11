import numpy as np

def demonstrate_fortress_problem():
    """
    Demonstrates that for n=4 guards on a sphere, there are unseen external points.
    """
    # 1. Define 4 guard positions on the unit sphere.
    # These are the vertices of a regular tetrahedron.
    p1 = np.array([0, 0, 1])
    p2 = np.array([2 * np.sqrt(2) / 3, 0, -1/3])
    p3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1/3])
    p4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1/3])

    guards = [p1, p2, p3, p4]
    print("Four guards are placed at the vertices of a regular tetrahedron:")
    for i, p in enumerate(guards):
        print(f"p{i+1}: {np.round(p, 3)}")
    print("-" * 30)

    # 2. Find a vertex of the polyhedron K defined by {x | p_i . x <= 1}.
    # A vertex is formed by the intersection of 3 tangent planes p_i . x = 1.
    # We solve the system of linear equations:
    # p1 . x = 1
    # p2 . x = 1
    # p3 . x = 1
    
    # Matrix A has guards as rows, b is the vector of 1s
    A = np.array([p1, p2, p3])
    b = np.array([1, 1, 1])

    # Solve Ax = b for the vertex x
    try:
        vertex = np.linalg.solve(A, b)
        print("Found a vertex of the polyhedron of unseen points.")
        print(f"Vertex v: {np.round(vertex, 3)}")
        print("-" * 30)
    except np.linalg.LinAlgError:
        print("The chosen guard vectors are not linearly independent.")
        return

    # 3. Check if this vertex is outside the unit sphere.
    norm = np.linalg.norm(vertex)
    print("Checking if the vertex is outside the unit ball (|v| > 1):")
    print(f"Magnitude |v|: {norm:.4f}")
    if norm > 1:
        print("Result: The vertex is indeed outside the unit ball.")
    else:
        print("Result: The vertex is not outside the unit ball. This configuration is not general.")
    print("-" * 30)

    # 4. Check if this vertex is seen by ANY of the guards.
    # It is unseen if p_i . v <= 1 for ALL guards.
    print("Checking if the vertex v is seen by any guard (p_i . v > 1?):")
    unseen = True
    for i, p in enumerate(guards):
        dot_product = np.dot(p, vertex)
        print(f"p{i+1} . v = {dot_product:.4f}")
        if dot_product > 1 + 1e-9: # Add tolerance for float precision
            unseen = False
    
    print("-" * 30)
    if unseen:
        print(f"Conclusion: The point v = {np.round(vertex, 3)} is outside the sphere and is not seen by any guard.")
        print("This demonstrates that 4 guards are not sufficient.")
        print("This logic extends to any finite number of guards.")
    else:
        print("Conclusion: The chosen vertex is seen by another guard, so it is not a blind spot.")


if __name__ == '__main__':
    demonstrate_fortress_problem()
