import sympy

def find_max_size_of_S():
    """
    This script finds the largest possible size for a set S of non-real eigenvalues
    of a matrix A satisfying A^3 = A^*.

    This is equivalent to finding the number of non-real solutions to the equation
    lambda^3 = conjugate(lambda).
    """

    # Define lambda as a complex number lambda = x + i*y
    x, y = sympy.symbols('x y', real=True)
    lam = x + sympy.I * y

    # Set up the equation lambda^3 = conjugate(lambda)
    equation = sympy.Eq(lam**3, sympy.conjugate(lam))

    # To solve this, we separate the equation into its real and imaginary parts.
    # The equation to solve is lambda**3 - conjugate(lambda) = 0
    real_part = sympy.re(equation.lhs - equation.rhs)
    imag_part = sympy.im(equation.lhs - equation.rhs)

    print("The condition on eigenvalues lambda leads to the equation: lambda**3 = conjugate(lambda).")
    print("Let lambda = x + iy. This gives a system of two polynomial equations for real x and y:")
    # The 'numbers in the final equation' are the coefficients and powers.
    print(f"1. Real part: {real_part} = 0")
    print(f"2. Imaginary part: {imag_part} = 0")

    # Solve the system of equations for (x, y)
    solutions = sympy.solve([real_part, imag_part], (x, y))

    # Filter for non-real solutions, which are those where y is not zero.
    non_real_solutions = []
    for sol_x, sol_y in solutions:
        if sol_y != 0:
            non_real_solutions.append(sol_x + sympy.I * sol_y)

    print("\nThe set of possible non-real eigenvalues is the set of non-real solutions:")
    print(non_real_solutions)
    
    # The largest size of S is the number of these possible non-real eigenvalues.
    max_size = len(non_real_solutions)

    print(f"\nA matrix A with these eigenvalues can be constructed (e.g., A = diag{tuple(non_real_solutions)}).")
    print(f"Thus, the largest possible size of the set S is {max_size}.")

if __name__ == '__main__':
    find_max_size_of_S()
    # The final answer is the integer value of the maximum size.
    # print(f"\n<<<2>>>")