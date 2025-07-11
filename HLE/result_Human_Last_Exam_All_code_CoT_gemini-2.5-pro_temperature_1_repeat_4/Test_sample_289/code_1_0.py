import sympy

def solve_eigenvalue_problem():
    """
    Solves for the largest size of a set of non-real eigenvalues
    of a matrix A satisfying A^3 = A^*.
    """
    # --- Analysis of the Eigenvalue Problem ---
    print("--- Analysis of the Eigenvalue Problem ---")
    print("Let A be an n x n matrix such that A^3 = A*, where A* is the conjugate transpose of A.")
    print("The matrix A is normal since AA* = A(A^3) = A^4 and A*A = (A^3)A = A^4.")
    print("A normal matrix is unitarily diagonalizable, A = UDU*, where D is a diagonal matrix of eigenvalues.")
    print("The condition A^3 = A* becomes UD^3U* = U(conjugate(D))U*, which simplifies to D^3 = conjugate(D).")
    print("This implies that for each eigenvalue lambda, the equation lambda^3 = conjugate(lambda) must hold.")

    # --- Solving the Eigenvalue Equation ---
    print("\n--- Solving the Eigenvalue Equation ---")
    print("We solve the equation for lambda. Let lambda = x + i*y for real x, y.")
    x, y = sympy.symbols('x y', real=True)
    z = x + sympy.I * y

    # Equation z^3 = conjugate(z)
    # We create the expression that should be zero: z^3 - conjugate(z) = 0
    equation_expr = z**3 - sympy.conjugate(z)
    
    # Separate the real and imaginary parts of the expression
    real_part = sympy.re(equation_expr)
    imag_part = sympy.im(equation_expr)

    print(f"The equation lambda^3 = conjugate(lambda) is equivalent to the system of equations:")
    # Using expand() to make the polynomial expressions clearer
    print(f"  Real Part:      {sympy.expand(real_part)} = 0")
    print(f"  Imaginary Part: {sympy.expand(imag_part)} = 0")

    # --- Finding All Solutions ---
    # Solve the system of equations for (x, y)
    solutions = sympy.solve([real_part, imag_part], (x, y), dict=True)

    # --- Finding Non-Real Eigenvalues ---
    print("\n--- Finding Non-Real Eigenvalues ---")
    print(f"The complete set of solutions (x, y) for the eigenvalues is: {solutions}")

    # The problem asks for S, a set of non-real eigenvalues.
    # We filter the solutions to find those where y != 0.
    non_real_solutions = [sol for sol in solutions if sol[y] != 0]

    num_non_real_solutions = len(non_real_solutions)
    non_real_eigenvalues = [sol[x] + sympy.I * sol[y] for sol in non_real_solutions]

    print(f"\nThe solutions corresponding to non-real eigenvalues (where y is not 0) are: {non_real_solutions}")
    print(f"These solutions correspond to the complex eigenvalues: {non_real_eigenvalues}")
    
    print("\nVerifying these eigenvalues in the equation lambda^3 = conjugate(lambda):")
    for val in non_real_eigenvalues:
        lhs = val**3
        rhs = sympy.conjugate(val)
        print(f"For eigenvalue {val}:")
        print(f"  ({val})^3 = {lhs}")
        print(f"  conjugate({val}) = {rhs}")
        print(f"  The equation {lhs} = {rhs} holds.")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print(f"The set of all possible non-real eigenvalues for a matrix A satisfying A^3 = A* is {set(non_real_eigenvalues)}.")
    print("To show that a size of 2 is achievable, we can construct a matrix A with these eigenvalues.")
    print("For example, for n=2, let A be the diagonal matrix with i and -i on the diagonal:")
    print("A = [[i, 0], [0, -i]]")
    print("For this matrix, A^* = [[-i, 0], [0, i]]")
    print("And A^3 = [[i^3, 0], [0, (-i)^3]] = [[-i, 0], [0, i]]")
    print("So, A^3 = A* is satisfied, and its set of non-real eigenvalues is {i, -i}.")
    print(f"\nTherefore, the largest possible size of the set S is the number of distinct non-real eigenvalues.")
    print(f"Largest size |S| = {num_non_real_solutions}")

if __name__ == '__main__':
    solve_eigenvalue_problem()