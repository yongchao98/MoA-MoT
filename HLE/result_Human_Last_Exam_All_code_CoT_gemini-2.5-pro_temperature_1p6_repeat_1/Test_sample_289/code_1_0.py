import sympy

def solve_eigenvalue_problem():
    """
    Solves for the eigenvalues lambda of a matrix A satisfying A^3 = A*
    and determines the maximum number of distinct non-real eigenvalues.
    """
    # The condition on the matrix A implies that each eigenvalue lambda must
    # satisfy the equation lambda^3 = conjugate(lambda).
    # Let lambda = x + i*y, where x and y are real. The equation becomes:
    # (x + i*y)^3 = x - i*y
    # This can be split into a system of two polynomial equations for x and y.
    
    x, y = sympy.symbols('x y', real=True)
    
    # Equation for the real parts
    # Real part of (x + iy)^3 is x^3 - 3xy^2
    eq1 = sympy.Eq(x**3 - 3*x*y**2, x)
    
    # Equation for the imaginary parts
    # Imaginary part of (x + iy)^3 is 3x^2y - y^3
    eq2 = sympy.Eq(3*x**2*y - y**3, -y)
    
    print("The problem is equivalent to finding the number of non-real solutions to lambda^3 = conjugate(lambda).")
    print("Let lambda = x + i*y. This gives the system of equations:")
    print(f"1)  {eq1.lhs} = {eq1.rhs}")
    print(f"2)  {eq2.lhs} = {eq2.rhs}")
    print("-" * 30)

    # Solve the system of equations for (x, y)
    solutions = sympy.solve([eq1, eq2], (x, y))
    
    print("The set of all possible eigenvalues (lambda = x + iy) is:")
    all_eigenvalues = []
    for sol in solutions:
        x_val, y_val = sol
        # Construct the complex number lambda
        lambda_val = x_val + sympy.I * y_val
        all_eigenvalues.append(lambda_val)
        print(f"-> {lambda_val}")
        
    # Filter for non-real eigenvalues (where the imaginary part y is non-zero)
    non_real_eigenvalues = [val for val in all_eigenvalues if sympy.im(val) != 0]
    
    print("\n" + "-" * 30)
    print("The subset S of non-real eigenvalues is:")
    for val in non_real_eigenvalues:
        print(f"-> {val}")

    # The size of this set is the largest possible size for S.
    max_size = len(non_real_eigenvalues)
    
    print("\n" + "-" * 30)
    print(f"The largest possible size of the set S of non-real eigenvalues is |S| = {max_size}.")

if __name__ == '__main__':
    solve_eigenvalue_problem()