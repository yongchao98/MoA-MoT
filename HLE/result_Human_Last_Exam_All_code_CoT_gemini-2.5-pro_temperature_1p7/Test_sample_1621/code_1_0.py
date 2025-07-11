import sympy

def solve_problem():
    """
    This function demonstrates the solution for n=4 by constructing the
    matrices based on quaternion algebra and calculating the determinant symbolically.
    """
    print("The problem is equivalent to finding for which n a real division algebra structure exists.")
    print("By the Bott-Milnor-Kervaire theorem, this is possible only for n = 1, 2, 4, and 8.")
    print("Thus, there are 4 such natural numbers.\n")
    print("--- Verification for n=4 (Quaternions) ---")

    # Define symbolic variables for the vector x
    x1, x2, x3, x4 = sympy.symbols('x1 x2 x3 x4')
    x = sympy.Matrix([x1, x2, x3, x4])

    # For n=4, we can construct the matrices A_i based on quaternion multiplication.
    # Let A_i be the matrix for right-multiplication by the i-th basis vector (1, i, j, k).
    # A_i * y represents y * e_i
    
    # A1 corresponds to right-multiplication by 1 (the identity)
    A1 = sympy.eye(4)
    
    # A2 corresponds to right-multiplication by i
    # (a+bi+cj+dk)i = -b+ai-dj+ck -> (-b, a, -d, c)
    A2 = sympy.Matrix([
        [0, -1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, -1],
        [0, 0, 1, 0]
    ])
    
    # A3 corresponds to right-multiplication by j
    # (a+bi+cj+dk)j = -c+di+aj-bk -> (-c, d, a, -b)
    A3 = sympy.Matrix([
        [0, 0, -1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, -1, 0, 0]
    ])

    # A4 corresponds to right-multiplication by k
    # (a+bi+cj+dk)k = -d-ci+bj+ak -> (-d, -c, b, a)
    A4 = sympy.Matrix([
        [0, 0, 0, -1],
        [0, 0, -1, 0],
        [0, 1, 0, 0],
        [1, 0, 0, 0]
    ])
    
    print("\nThe matrices A_i for n=4 are:")
    print("A1:"); sympy.pprint(A1)
    print("A2:"); sympy.pprint(A2)
    print("A3:"); sympy.pprint(A3)
    print("A4:"); sympy.pprint(A4)

    # Calculate the column vectors A_i * x
    c1 = A1 * x
    c2 = A2 * x
    c3 = A3 * x
    c4 = A4 * x

    # Construct the matrix B(x) = [A1*x | A2*x | A3*x | A4*x]
    B = sympy.Matrix.hstack(c1, c2, c3, c4)

    print("\nThe matrix B(x) is:")
    sympy.pprint(B)

    # Calculate the determinant of B(x)
    det_B = B.det()
    
    # Simplify the expression for the determinant
    det_B_simplified = sympy.simplify(det_B)

    print("\nThe determinant of B(x) is calculated as:")
    # Using sympy.pretty to format the equation nicely
    final_equation = sympy.Eq(sympy.Symbol('det(B(x))'), det_B_simplified)
    sympy.pprint(final_equation)

    print("\nThis determinant is zero if and only if x1=x2=x3=x4=0 (i.e., x is the zero vector).")
    print("For any non-zero x, the determinant is non-zero, so the vectors are linearly independent.")

solve_problem()
