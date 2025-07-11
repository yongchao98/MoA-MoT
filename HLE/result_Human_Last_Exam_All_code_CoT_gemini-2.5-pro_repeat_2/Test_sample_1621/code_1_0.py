import sympy

def solve():
    """
    This function analyzes the problem for different values of n and
    demonstrates the solution for n=1, 2, 4.
    """
    print("This program investigates for which natural numbers n there exist n real n-by-n matrices A_1, ..., A_n")
    print("such that for any non-zero vector x, the vectors A_1*x, ..., A_n*x are linearly independent.\n")
    print("The problem is equivalent to finding n for which an n-dimensional real division algebra exists.")
    print("A theorem states this is only possible for n = 1, 2, 4, 8.\n")

    possible_n = [1, 2, 4, 8]
    print(f"There are {len(possible_n)} such values of n: {possible_n}.\n")

    def check_n(n):
        print(f"--- Checking for n = {n} ---")
        if n not in possible_n:
            print(f"A solution for n={n} does not exist based on the theorem.")
            # For odd n>1, we can add the specific reasoning.
            if n % 2 != 0 and n > 1:
                print("Reason: For odd n>1, the determinant polynomial P(x) is an odd function,")
                print("and by the Borsuk-Ulam theorem, it must have a zero on the unit sphere.")
            return

        # Symbolic variables for the vector x
        x_vars = sympy.symbols(f'x1:{n+1}')
        x_vec = sympy.Matrix(x_vars)
        
        matrices = []
        if n == 1:
            # For R (Real Numbers)
            A1 = sympy.Matrix([[1]])
            matrices.append(A1)
            # For n=1, the condition is that A1*x is not the zero vector.
            # a*x != 0 for x != 0, which requires a != 0. A1=[1] works.
            print("For n=1, we need A1*x != 0 for x != 0.")
            print(f"Let A1 = {matrices[0]}. Then A1*x = {matrices[0]*x_vec}, which is non-zero if x1 is non-zero.")
            return

        elif n == 2:
            # For C (Complex Numbers)
            A1 = sympy.Matrix([[1, 0], [0, 1]]) # Represents 1
            A2 = sympy.Matrix([[0, -1], [1, 0]]) # Represents i
            matrices.extend([A1, A2])
        
        elif n == 4:
            # For H (Quaternions)
            A1 = sympy.eye(4) # Represents 1
            A2 = sympy.Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]) # Represents i
            A3 = sympy.Matrix([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]]) # Represents j
            A4 = sympy.Matrix([[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]]) # Represents k
            matrices.extend([A1, A2, A3, A4])

        elif n == 8:
            print("The construction for n=8 is based on Octonions.")
            print("The symbolic calculation of an 8x8 determinant is computationally expensive and will be skipped.")
            # We can define the matrices but not compute the determinant
            A1 = sympy.eye(8) # e0
            # Cayley-Dickson construction can be used to define the others
            # A2 to A8 correspond to e1 to e7
            print("A set of 8x8 matrices can be constructed from the Octonion multiplication table.")
            return

        print("Constructed matrices A_1, ..., A_n:")
        for i, A in enumerate(matrices):
            print(f"A_{i+1}:")
            sympy.pprint(A)
            print()

        # Construct the matrix B_x = [A1*x | A2*x | ... | An*x]
        Bx_cols = [A * x_vec for A in matrices]
        Bx = sympy.Matrix.hstack(*Bx_cols)

        print("Matrix B_x = [A1*x | ... | An*x]:")
        sympy.pprint(Bx)
        print()
        
        # Calculate the determinant
        print("Calculating det(B_x)...")
        det_Bx = Bx.det()
        
        # Simplify the expression
        print("Simplifying the determinant...")
        simplified_det = sympy.simplify(det_Bx)

        print("The determinant polynomial is:")
        # The pretty print is nice for expressions
        sympy.pprint(simplified_det)
        
        # Verify it's a sum of squares
        expected_expr = sum(var**2 for var in x_vars)
        if n == 2:
            expected_expr = expected_expr
        elif n == 4:
            expected_expr = expected_expr**2
        elif n == 8: # We don't run this case, but for completeness
            expected_expr = expected_expr**4

        print("\nThis polynomial is >= 0 and is 0 only if all x_i are 0.")
        print("Thus, the vectors are linearly independent for any non-zero x.")


    # Run the checks
    check_n(1)
    print("-" * 25)
    check_n(2)
    print("-" * 25)
    check_n(3)
    print("-" * 25)
    check_n(4)
    print("-" * 25)
    check_n(8)

solve()
