import sympy as sp

def solve_normal_cone():
    """
    This function calculates the explicit representation of the normal cone for the given problem.
    """
    # Step 0: Define variables and functions
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    x = sp.Matrix([x1, x2, x3])
    x_star_val = sp.Matrix([2, 0, -1])

    # Constraint functions g_i(x)
    g_exprs = [
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ]

    print("Feasible set F is defined by g_i(x) <= 0 for i=1,2,3,4.")
    print("The point is x* = (2, 0, -1).\n")

    # Step 1: Identify active constraints
    print("Step 1: Identify active constraints by evaluating g_i(x*)")
    active_indices = []
    active_constraints_info = []
    for i, g_expr in enumerate(g_exprs):
        val = g_expr.subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
        status = "active" if sp.Abs(val) < 1e-9 else "inactive"
        if status == "active":
            active_indices.append(i)
            active_constraints_info.append(f"g_{i+1}(x*) = {val} (active)")
        else:
            active_constraints_info.append(f"g_{i+1}(x*) = {val} (inactive)")
    
    print("\n".join(active_constraints_info))
    print(f"\nThe active constraints are g_i for i in {[i + 1 for i in active_indices]}.\n")

    # Step 2: Check for convexity (conceptually, not coded, but stated)
    print("Step 2: Check for convexity")
    print("The constraint functions g_1, g_2 are convex (quadratic with positive semi-definite Hessians).")
    print("The constraint functions g_3, g_4 are linear, hence convex.")
    print("Since all constraint functions are convex, the feasible set F is convex.")
    print("For convex sets, the normal cone is the conic hull of the gradients of active constraints.\n")

    # Step 3: Compute gradients of active constraints
    print("Step 3: Compute gradients of active constraints at x*")
    active_gradients = []
    for i in active_indices:
        g_expr = g_exprs[i]
        grad_g = sp.Matrix([sp.diff(g_expr, var) for var in x])
        grad_g_at_x_star = grad_g.subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
        active_gradients.append(grad_g_at_x_star)
        print(f"Gradient of g_{i+1} at x*: {list(grad_g_at_x_star)}")
    print("")

    # Step 4: Formulate the normal cone
    print("Step 4: Formulate the normal cone T_F^°(x*)")
    print("A vector s = (s1, s2, s3) is in the normal cone if it's a non-negative linear combination of the active gradients.")
    
    grad_vectors_str = [f"lambda_{i+1}*{list(active_gradients[i])}" for i in range(len(active_gradients))]
    print(f"s = {' + '.join(grad_vectors_str)}, with lambda_1, lambda_2, lambda_3 >= 0.")
    print("\nLet's write out the components of s:")
    
    # Printing the equation with each number
    s1_eq = f"s1 = lambda_1*({active_gradients[0][0]}) + lambda_2*({active_gradients[1][0]}) + lambda_3*({active_gradients[2][0]})"
    s2_eq = f"s2 = lambda_1*({active_gradients[0][1]}) + lambda_2*({active_gradients[1][1]}) + lambda_3*({active_gradients[2][1]})"
    s3_eq = f"s3 = lambda_1*({active_gradients[0][2]}) + lambda_2*({active_gradients[1][2]}) + lambda_3*({active_gradients[2][2]})"
    print(s1_eq)
    print(s2_eq)
    print(s3_eq)
    print("")

    # Step 5: Simplify and represent
    print("Step 5: Simplify to find the explicit representation")
    print("Simplifying the equations for the components of s:")
    print("s1 = 2*lambda_1 - 2*lambda_2")
    print("s2 = 0")
    print("s3 = lambda_3")

    print("\nFrom these equations, we deduce the conditions on s1, s2, s3:")
    print("1. For s1 = 2*(lambda_1 - lambda_2), since lambda_1, lambda_2 >= 0, s1 can take any real value.")
    print("2. The component s2 is always 0.")
    print("3. For s3 = lambda_3, since lambda_3 >= 0, we must have s3 >= 0.")
    print("\nTherefore, the explicit representation of the normal cone is:")
    print("T_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s2 = 0, s3 >= 0 }")

if __name__ == '__main__':
    solve_normal_cone()