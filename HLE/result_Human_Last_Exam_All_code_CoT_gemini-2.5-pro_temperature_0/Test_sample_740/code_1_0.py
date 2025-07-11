import numpy as np
import sys

def solve_for_n(n):
    """
    Solves the integral equation for a given integer n.
    Finds the values of lambda for which there is no solution.
    """
    if not isinstance(n, int) or n < 1:
        print("n must be a positive integer.", file=sys.stderr)
        return

    print(f"--- Results for n = {n} ---")

    # The problem reduces to analyzing the eigenvalues of the matrix A,
    # where A[i, j] = 1 / (n + j - i) for i, j in {0, ..., n-1}.
    # This matrix represents the integral operator in the basis {1, x, ..., x^(n-1)}.
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = 1.0 / (n + j - i)

    # The characteristic values lambda are the roots of the polynomial
    # det(I - lambda*A) = 0, which is equivalent to det(A - (1/lambda)*I) = 0.
    # Let mu = 1/lambda. The characteristic polynomial for mu is det(A - mu*I) = 0.
    # Let p(mu) = c_n*mu^n + ... + c_0.
    # The polynomial for lambda is c_0*lambda^n + ... + c_n = 0.
    poly_coeffs_mu = np.poly(A)
    poly_coeffs_lambda = poly_coeffs_mu[::-1]

    # Print the equation for lambda
    equation_str = []
    for i, coeff in enumerate(poly_coeffs_lambda):
        power = n - i
        if abs(coeff) > 1e-9:
            sign = "-" if coeff < 0 else "+"
            if i == 0:
                sign = "" if coeff > 0 else "-"
            
            coeff_val = abs(coeff)
            
            if power > 1:
                term = f"{sign} {coeff_val:.4f}*λ^{power}"
            elif power == 1:
                term = f"{sign} {coeff_val:.4f}*λ"
            else:
                term = f"{sign} {coeff_val:.4f}"
            equation_str.append(term)
    
    print("The equation for λ is:")
    print(" ".join(equation_str).lstrip('+ ').lstrip(), "= 0")


    # According to the Fredholm alternative, if lambda = 1/mu for an eigenvalue mu,
    # a solution exists if and only if the function f(x)=1 is orthogonal to the
    # eigenfunctions of the adjoint operator corresponding to mu.
    # Since the kernel is symmetric, the operator is self-adjoint.
    # The condition for a solution to exist is integral(1 * phi(x) dx) = 0,
    # where phi is the eigenfunction for mu.
    
    # We find the values of lambda for which there is NO solution.
    # This occurs when the orthogonality condition is NOT met.

    # Eigenvalues mu and right eigenvectors p of A.
    # The columns of p are the eigenvectors.
    mu_values, p_vectors = np.linalg.eig(A)

    # The vector b comes from the integral of the basis functions used for phi.
    # phi(x) = p_0*x^0 + p_1*x^1 + ...
    # integral(phi(x) dx) = p_0/1 + p_1/2 + ... = b.dot(p)
    b = np.array([1.0 / (i + 1) for i in range(n)])

    print("\nValues of λ for which no solution exists:")
    
    no_solution_lambdas = []
    for i in range(n):
        mu = mu_values[i]
        p = p_vectors[:, i]
        
        # Check orthogonality condition: integral(1 * phi(x) dx) != 0
        dot_product = np.dot(b, p)
        
        # We check if the dot product is significantly different from zero.
        if abs(dot_product) > 1e-9:
            lambda_val = 1.0 / mu
            no_solution_lambdas.append(lambda_val)
            
    # Sort the values for consistent output
    no_solution_lambdas.sort(key=lambda x: (x.real, x.imag))
    for val in no_solution_lambdas:
        if abs(val.imag) < 1e-9:
            print(f"{val.real:.8f}")
        else:
            print(f"{val.real:.8f} + {val.imag:.8f}i")
    print("-" * 25 + "\n")


if __name__ == '__main__':
    # We solve the problem for n = 1, 2, and 3 as examples.
    solve_for_n(1)
    solve_for_n(2)
    solve_for_n(3)