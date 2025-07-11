import numpy as np

# Set the parameter n for the equation.
# For n=1, the only value is λ=1.
# We choose n=2 for a non-trivial example.
n = 2

# Based on Fredholm theory, the problem reduces to an eigenvalue problem.
# The values of λ for which there are no solutions are the reciprocals
# of the eigenvalues of the matrix A, defined below.

# The matrix A is defined as A_jk = 1 / (n - j + k)
# where j, k are 0-indexed, from 0 to n-1.
A = np.zeros((n, n))
for j in range(n):
    for k in range(n):
        A[j, k] = 1.0 / (n - j + k)

# The characteristic equation for λ is derived from det(A - (1/λ)I) = 0.
# This results in a polynomial in λ. The coefficients of this polynomial
# are the reverse of the coefficients of the characteristic polynomial of A.
# np.poly(A) finds coefficients of det(xI - A)=0. Let x=mu.
mu_poly_coeffs = np.poly(A)
lambda_poly_coeffs = mu_poly_coeffs[::-1]

# The values of λ are the roots of this polynomial. They can also be
# calculated directly as the reciprocals of the eigenvalues of A.
mu_values = np.linalg.eigvals(A)
lambda_values = 1.0 / mu_values

# Print the results, including the equation and its constituent numbers.
print(f"For n={n}, the characteristic equation for λ, whose roots are the values where no solution exists, is:")

# Build and print the equation string from its coefficients
equation_parts = []
n_poly = len(lambda_poly_coeffs) - 1 # Degree of the polynomial
for i, coeff in enumerate(lambda_poly_coeffs):
    power = n_poly - i
    
    if abs(coeff) < 1e-9:
        continue # Skip terms with near-zero coefficients
    
    coeff_str = f"{abs(coeff):.4f}"
    
    if power > 1:
        term = f"{coeff_str}*λ^{power}"
    elif power == 1:
        term = f"{coeff_str}*λ"
    else: # power == 0
        term = f"{coeff_str}"
        
    # Prepend the sign for the term
    # The first term's sign is handled slightly differently
    if i > 0 and equation_parts:
         equation_parts.append(f" {'-' if coeff < 0 else '+'} ")
    elif coeff < 0:
         equation_parts.append("-")
    equation_parts.append(term)

final_equation_str = "".join(equation_parts) + " = 0"
print(final_equation_str)

print("\nThe numbers in the final equation (coefficients from highest power to constant term):")
for c in lambda_poly_coeffs:
    print(f"{c:.4f}")
    
print("\nThe values of λ (roots of the equation):")
# The eigenvalues of the underlying operator are real. Numerical results may have small imaginary parts.
# We sort the real parts for consistent output.
for val in sorted(lambda_values.real):
    print(f"{val:.4f}")
