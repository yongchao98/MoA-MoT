import numpy as np

# We need to find the average of the complex coordinates z which are the singularities of the system.
# These singularities are the roots of the polynomial in the denominator of the B-field equation.
# The polynomial is P(z) = 4*z^4 - z^3 + z^2 + 1 = 0.
# The coefficients are for z^4, z^3, z^2, z^1, z^0 respectively.
coeffs = [4, -1, 1, 0, 1]

# The degree of the polynomial is the number of roots.
n = 4

# Using Vieta's formulas, the sum of the roots is -a_{n-1}/a_n.
# Here, a_n is the coefficient of z^4 (coeffs[0]) and a_{n-1} is the coefficient of z^3 (coeffs[1]).
sum_of_roots = -coeffs[1] / coeffs[0]

# The average of the roots is the sum divided by the number of roots.
average_of_roots = sum_of_roots / n

# The final equation is Average = Sum / Number of Roots.
# The prompt requests to output each number in this equation.
print("The final equation for the average of the coordinates is:")
print(f"Average = (Sum of Roots) / (Number of Roots)")
print(f"Average = ({sum_of_roots}) / ({n})")
print(f"Average = {average_of_roots}")

# As a verification, we can compute the roots numerically.
roots = np.roots(coeffs)
print("\n--- Numerical Verification ---")
print(f"The 4 complex roots are:")
for root in roots:
    print(root)
print(f"Numerically computed average: {np.mean(roots).real}")