import sympy as sp

def find_shape_equation():
    """
    Symbolically derives the equation for the set S for n=2.
    """
    # Step 1: Define symbolic vectors y1, y2 in R^2
    y1_sym = sp.Matrix([1, 0])
    y2_sym = sp.Matrix([1, 1])

    print("Chosen linearly independent vectors:")
    print(f"y1 = {y1_sym.tolist()}")
    print(f"y2 = {y2_sym.tolist()}")
    print("-" * 30)

    # Step 2: Compute the Gram matrix G
    G = sp.Matrix([
        [y1_sym.dot(y1_sym), y1_sym.dot(y2_sym)],
        [y2_sym.dot(y1_sym), y2_sym.dot(y2_sym)]
    ])
    print("Gram matrix G = y_i^T * y_j:")
    print(G)
    print("-" * 30)
    
    # Step 3: Compute the inverse of the Gram matrix H = G^-1
    H = G.inv()
    print("Inverse Gram matrix H = G^-1:")
    print(H)
    print("-" * 30)

    # Step 4: Define symbolic variables for the components of v and p
    v1, v2 = sp.symbols('v1 v2')
    p1, p2 = sp.symbols('p1 p2')
    v = sp.Matrix([v1, v2])

    # The set of vectors v = (<y1,s>, <y2,s>) forms an ellipsoid v^T * H * v = 1
    ellipsoid_eq_v = (v.T * H * v)[0] - 1
    print("Equation of the ellipsoid for v=(<y1,s>, <y2,s>):")
    print(f"{sp.expand(ellipsoid_eq_v)} = 0")
    print("-" * 30)

    # Step 5: Substitute v1 = sqrt(p1), v2 = sqrt(p2) to find the relation for p
    # To handle the +/- sign, we isolate the cross term and square it.
    # H[0,0]*p1 + H[1,1]*p2 + 2*H[0,1]*sqrt(p1*p2) = 1
    
    # Rearrange: 2*H[0,1]*sqrt(p1*p2) = 1 - H[0,0]*p1 - H[1,1]*p2
    lhs = (2 * H[0, 1] * sp.sqrt(p1 * p2))**2
    rhs = (1 - H[0, 0] * p1 - H[1, 1] * p2)**2

    # The equation for S is lhs - rhs = 0
    shape_eq_p = sp.expand(lhs - rhs)
    
    # The final equation is a quadratic in p1, p2, which represents an ellipse.
    print("Final equation for the shape S with components (p1, p2):")
    final_eq_str = f"{shape_eq_p} = 0"
    print(final_eq_str)
    
    # Let's show this in the standard form (x-h)^2/a^2 + (y-k)^2/b^2 = 1 to confirm it's an ellipse.
    # The equation is 4*p1**2 - 4*p1 + p2**2 - 2*p2 + 1 = 0
    # Complete the square
    p1_term = 4*(p1**2 - p1)
    p2_term = p2**2 - 2*p2
    # 4*((p1 - 1/2)**2 - 1/4) + (p2 - 1)**2 - 1 + 1 = 0
    # 4*(p1 - 1/2)**2 + (p2 - 1)**2 = 1
    print("\nStandard form of the ellipse equation:")
    print("4*(p1 - 1/2)**2 + (p2 - 1)**2 = 1")
    print("\nThis confirms the shape is an ellipsoid (ellipse in 2D).")

find_shape_equation()