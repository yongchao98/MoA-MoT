import sympy

def solve():
    """
    This function demonstrates the existence of the required matrices for n=1, 2, and 4
    by constructing them and computing the determinant polynomial P(x).
    The construction is based on real division algebras (Reals, Complex, Quaternions).
    """
    possible_n = [1, 2, 4]
    
    print("This script checks for which n, there exist n matrices A_i such that")
    print("det([A_1*x, ..., A_n*x]) is non-zero for any non-zero vector x.\n")

    # Case n=1 (Real numbers)
    n = 1
    x_vars = sympy.symbols(f'x_1:{n+1}')
    x = sympy.Matrix(x_vars)
    A1 = sympy.Matrix([[1]])
    matrices = [A1]
    Bx = sympy.Matrix.hstack(*(A * x for A in matrices))
    P_x = Bx.det()
    print(f"For n = {n}:")
    print(f"The determinant polynomial is: {P_x}")
    print("This is non-zero for any non-zero x.\n")
    
    # Case n=2 (Complex numbers)
    n = 2
    x_vars = sympy.symbols(f'x_1:{n+1}')
    x = sympy.Matrix(x_vars)
    # A1 corresponds to multiplication by 1
    A1 = sympy.eye(n)
    # A2 corresponds to multiplication by i
    A2 = sympy.Matrix([[0, -1], [1, 0]])
    matrices = [A1, A2]
    Bx = sympy.Matrix.hstack(*(A * x for A in matrices))
    P_x = sympy.simplify(Bx.det())
    print(f"For n = {n}:")
    print(f"The determinant polynomial is: {P_x}")
    print("This is a sum of squares, non-zero for any non-zero x.\n")
    
    # Case n=4 (Quaternions)
    n = 4
    x_vars = sympy.symbols(f'x_1:{n+1}')
    x = sympy.Matrix(x_vars)
    # A1 corresponds to multiplication by 1
    A1 = sympy.eye(n)
    # A2 corresponds to multiplication by i
    A2 = sympy.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    # A3 corresponds to multiplication by j
    A3 = sympy.Matrix([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0]])
    # A4 corresponds to multiplication by k
    A4 = sympy.Matrix([[0, 0, 0, -1], [0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
    matrices = [A1, A2, A3, A4]
    Bx = sympy.Matrix.hstack(*(A * x for A in matrices))
    P_x = sympy.simplify(Bx.det())
    print(f"For n = {n}:")
    print(f"The determinant polynomial is: {P_x}")
    print("This is a sum of squares squared, non-zero for any non-zero x.\n")
    
    print("The argument based on the Intermediate Value Theorem rules out all odd n > 1.")
    print("A deep theorem from algebraic topology shows that n can only be 1, 2, 4, or 8.")
    
    final_count = 4
    print(f"Thus, there are {final_count} such natural numbers (1, 2, 4, 8).")

solve()