import numpy as np
import sys

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the given integral equation has no solutions
    for a specified integer n >= 1.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be an integer greater than or equal to 1.", file=sys.stderr)
        return

    # For n=1, the equation u(x) = 1 + λ ∫u(y)dy leads to λ=1 having no solution.
    if n == 1:
        print("For n = 1, the equation is λ^1 - 1 = 0")
        print("The value of λ for which there is no solution is:")
        print("λ = 1.0")
        return

    # The problem can be converted to a matrix equation (I - λA)c = g.
    # No solution exists for λ = 1/μ where μ is an eigenvalue of A,
    # provided a condition on the forcing term g holds (which it does).
    # The matrix A has elements A_ij = 1/(n+j-i) for i, j from 0 to n-1.
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = 1.0 / (n + j - i)

    # The values of λ are the roots of det(I - λA) = 0.
    # Let p(μ) = det(μI - A) be the characteristic polynomial of A. The roots are eigenvalues μ.
    # We have det(I - λA) = λ^n * det((1/λ)I - A) = λ^n * p(1/λ).
    # If p(μ) = p_n*μ^n + p_{n-1}*μ^{n-1} + ... + p_0 (with p_n=1),
    # then λ^n * p(1/λ) = p_n + p_{n-1}λ + ... + p_0*λ^n.
    # The coefficients of the polynomial in λ are the reversed coefficients of the polynomial in μ.
    
    # Get coefficients of the characteristic polynomial of A, p(μ) = det(μI - A)
    # The coefficients are ordered from highest power to lowest, i.e., [p_n, p_{n-1}, ..., p_0]
    poly_coeffs_mu = np.poly(A)

    # Reverse to get coefficients of the polynomial in λ.
    # These coefficients correspond to powers 0, 1, ..., n.
    poly_coeffs_lambda = poly_coeffs_mu[::-1]

    print(f"For n = {n}, the values of λ for which the equation has no solution are the roots of the following polynomial equation:")

    equation_parts = []
    # Iterate through coefficients for powers λ^0, λ^1, ... λ^n
    for i, coeff in enumerate(poly_coeffs_lambda):
        power = i
        if np.isclose(coeff, 0):
            continue
        
        # Determine the sign and value string
        sign = ""
        if len(equation_parts) > 0:
            sign = "+ " if coeff > 0 else "- "
        elif coeff < 0:
            sign = "-"
        
        val_str = f"{abs(coeff):.4f}"
        if np.isclose(abs(coeff), 1.0):
            val_str = "" # Don't print 1.0 if it's a coefficient

        # Format variable part
        if power == 0:
            var_str = ""
            if val_str == "": # case for coeff 1 or -1 with no var
                val_str = "1.0000"
        elif power == 1:
            var_str = "λ"
        else:
            var_str = f"λ^{power}"

        equation_parts.append(f"{sign}{val_str}{var_str}")

    # Join parts and make it pretty
    final_equation = " ".join(equation_parts).strip()
    # Correcting leading "+ " if it exists
    if final_equation.startswith('+ '):
        final_equation = final_equation[2:]

    print(f"{final_equation} = 0")
    
    # The numerical solutions for λ are the reciprocals of the eigenvalues of A.
    eigenvalues_mu = np.linalg.eigvals(A)
    lambdas = 1.0 / eigenvalues_mu
    
    print("\nThe approximate numerical values for λ are:")
    # Sort for consistent output order
    for l in sorted(lambdas.real):
        print(f"λ = {l:.8f}")

# The prompt gives a general 'n', so we run the solution for a representative small case, n=2.
# You can change this value to solve for other n.
n_value = 2
solve_for_lambda(n_value)