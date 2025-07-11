import numpy as np

def solve_fortress_problem():
    """
    Solves the fortress problem for a unit ball in 3D by demonstrating
    that 4 guards are sufficient.
    """
    print("This program demonstrates that 4 guards are sufficient to watch the exterior of a sphere.")
    print("The problem is equivalent to finding the minimum number of points on a sphere")
    print("whose convex hull contains the sphere's center (the origin).\n")

    print("We argue that 3 guards are insufficient, as their convex hull is a triangle (a 2D object), which has no 3D interior.")
    print("We now show that 4 guards are sufficient by placing them at the vertices of a regular tetrahedron centered at the origin.\n")

    # 1. Define the vertices of a regular tetrahedron inscribed in the unit sphere.
    # These points will be our guard positions.
    g1 = np.array([0.0, 0.0, 1.0])
    g2 = np.array([2 * np.sqrt(2) / 3, 0.0, -1.0/3.0])
    g3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1.0/3.0])
    g4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1.0/3.0])

    guards = [g1, g2, g3, g4]
    
    print("Guard positions (vertices of a regular tetrahedron on the unit sphere):")
    for i, g in enumerate(guards):
        # Format the output for better readability
        print(f"  Guard {i+1}: ({g[0]: .4f}, {g[1]: .4f}, {g[2]: .4f})")

    # 2. Show that the origin is the centroid of these points.
    # For a regular tetrahedron centered at the origin, the origin is the centroid,
    # meaning it's the average of the vertex positions.
    # Origin = (1/4)*g1 + (1/4)*g2 + (1/4)*g3 + (1/4)*g4
    
    coefficient = 0.25
    combination_result = coefficient * (g1 + g2 + g3 + g4)
    
    print("\nTo show the origin is inside their convex hull, we express it as a convex combination:")
    print("Origin = c*g1 + c*g2 + c*g3 + c*g4, with c = 1/4.")
    print("\nCalculation:")
    print(f"{coefficient} * ({g1[0]:.4f}, {g1[1]:.4f}, {g1[2]:.4f}) +")
    print(f"{coefficient} * ({g2[0]:.4f}, {g2[1]:.4f}, {g2[2]:.4f}) +")
    print(f"{coefficient} * ({g3[0]:.4f}, {g3[1]:.4f}, {g3[2]:.4f}) +")
    print(f"{coefficient} * ({g4[0]:.4f}, {g4[1]:.4f}, {g4[2]:.4f})")
    
    print(f"\n= ({combination_result[0]:.4f}, {combination_result[1]:.4f}, {combination_result[2]:.4f})")

    print("\nThe result is the zero vector (the origin), proving the origin is inside the convex hull of the guard positions.")
    print("Since 3 guards are insufficient and 4 are sufficient, the minimum number is 4.")

solve_fortress_problem()