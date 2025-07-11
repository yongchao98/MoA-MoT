import sympy

def solve_for_lambda(n):
    """
    For a given integer n, this function calculates the values of λ for which
    the integral equation u(x) = 1 + λ ∫[0,1] dy (xⁿ - yⁿ)/(x - y) u(y)
    has no solutions.
    """
    try:
        n = int(n)
        if n < 1:
            raise ValueError
    except (ValueError, TypeError):
        print("Error: Please enter a positive integer for n.")
        return

    print(f"\nAnalyzing the integral equation for n = {n}:")
    
    # The problem reduces to finding λ = 1/μ, where μ are the eigenvalues of an n x n matrix B.
    # B[j, l] = 1 / (n + l - j) for j, l in {0, ..., n-1}.
    # We use sympy for precision as this matrix can be ill-conditioned.

    mu = sympy.Symbol('μ')
    
    # Construct the matrix B using rational numbers
    B = sympy.Matrix(n, n, lambda j, l: sympy.Rational(1, n + l - j))

    # The eigenvalues μ are the roots of the characteristic polynomial det(B - μ*I) = 0
    try:
        char_poly = B.charpoly(mu)
        # This part fulfills the requirement: "output each number in the final equation!"
        # by printing the characteristic equation for μ.
        equation_str = str(sympy.expand(char_poly.as_expr())).replace('**', '^')
        print("\nThe values of λ are the reciprocals of μ, which are the roots of the characteristic equation:")
        print(f"{equation_str} = 0\n")
    except Exception as e:
        print(f"Could not compute the characteristic polynomial: {e}")

    # Find eigenvalues (μ) and eigenvectors (d) of B
    try:
        eigen_data = B.eigenvects()
    except Exception as e:
        print(f"Could not compute eigenvalues. The matrix for n={n} may be too complex. Error: {e}")
        return
    
    # Construct vector b for the solvability check (b_j = 1/(j+1))
    b = sympy.Matrix(n, 1, lambda j, _: sympy.Rational(1, j + 1))
    
    lambdas_no_solution = []
    
    # For each eigenvalue, check the Fredholm alternative condition.
    # No solution exists if d.T * b != 0, where d is an eigenvector.
    for mu_val, multiplicity, e_vectors in eigen_data:
        # We only need to check one eigenvector for each eigenvalue
        d = e_vectors[0]
        dot_product = (d.T * b)[0]
        
        # Check if the dot product is non-zero
        if dot_product.simplify() != 0:
            # If non-zero, no solution exists for λ = 1/μ
            for _ in range(multiplicity):
                 lambdas_no_solution.append(1 / mu_val)

    if not lambdas_no_solution:
        print("No values of λ were found for which the equation has NO solution.")
        print("This may mean that for any λ where a unique solution does not exist, there are infinitely many solutions.")
    else:
        print("The values of λ for which the equation has no solutions are:")
        for lam in lambdas_no_solution:
            # Print the value in both exact (symbolic) and numerical form
            print(f"  λ = {sympy.pretty(lam, use_unicode=True)}")
            print(f"    (Numerical value: λ ≈ {lam.evalf(12)})\n")

if __name__ == '__main__':
    try:
        n_input = input("Enter the integer value for n (e.g., 2): ")
        solve_for_lambda(n_input)
    except KeyboardInterrupt:
        print("\nExecution aborted.")
