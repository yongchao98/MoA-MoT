import sympy

def solve_eigenvalue_problem():
    """
    This function solves for the largest possible size of a set of non-real
    eigenvalues of a matrix A satisfying A^3 = A^*.
    """
    print("The analysis shows that any eigenvalue lambda must satisfy the equation: lambda^3 = conjugate(lambda).")
    print("We will solve this equation to find all possible eigenvalues.\n")

    # Step 1: Set up the equation for eigenvalues.
    # Let lambda = x + i*y, where x and y are real.
    print("Step 1: Set up the equation for lambda = x + iy")
    x, y = sympy.symbols('x y', real=True)
    lambda_expr = x + sympy.I * y
    
    equation = sympy.Eq(lambda_expr**3, sympy.conjugate(lambda_expr))
    print(f"The equation is: (x + I*y)^3 = x - I*y\n")

    # Step 2: Separate the equation into its real and imaginary parts.
    print("Step 2: Separate into real and imaginary part equations")
    re_eq = sympy.Eq(sympy.re(equation.lhs), sympy.re(equation.rhs))
    im_eq = sympy.Eq(sympy.im(equation.lhs), sympy.im(equation.rhs))
    
    # Sympy can expand these directly.
    real_part_eq = sympy.expand_complex(re_eq)
    imag_part_eq = sympy.expand_complex(im_eq)
    
    print(f"Real part:      {real_part_eq}")
    print(f"Imaginary part: {imag_part_eq}\n")

    # Step 3: Solve the system of two equations for (x, y).
    print("Step 3: Solve the system of equations for (x, y)")
    solutions = sympy.solve([real_part_eq, imag_part_eq], (x, y))
    print(f"The real solutions for (x, y) are: {solutions}\n")

    # Step 4: Convert the solutions back to complex eigenvalues.
    print("Step 4: Determine all possible complex eigenvalues")
    all_eigenvalues = {sol[0] + sympy.I * sol[1] for sol in solutions}
    print(f"The set of all possible eigenvalues is: {all_eigenvalues}\n")

    # Step 5: Filter for the non-real eigenvalues.
    print("Step 5: Identify the non-real eigenvalues")
    non_real_eigenvalues = {val for val in all_eigenvalues if sympy.im(val) != 0}
    print("The set S of all possible non-real eigenvalues is the subset of the above with non-zero imaginary parts.")
    print(f"S = {non_real_eigenvalues}")
    print("\nThe members of this set are:")
    for val in non_real_eigenvalues:
        # In this final equation, we print each number as requested
        print(f"lambda = {val}")
    print()

    # Step 6: The size of this set is the maximum size |S|.
    # We can construct a matrix A = diag(i, -i) that fulfills the property
    # and has these two non-real eigenvalues.
    max_size = len(non_real_eigenvalues)
    print("Step 6: Conclusion")
    print("The largest possible size |S| is the number of unique non-real eigenvalues found.")
    print(f"The largest size |S| is {max_size}.")

if __name__ == '__main__':
    solve_eigenvalue_problem()