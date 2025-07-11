# Based on the analysis above, we have derived the polynomials for the difference.
# The difference is p(d)/q(d).
# p(d) = -3d^2 + 3d
# q(d) = -2d^2 + 2d + 4

# Let's define the coefficients of these polynomials.
# For p(d), the coefficients for d^2, d^1, and d^0 are -3, 3, and 0.
p_coeffs = [-3, 3, 0]

# For q(d), the coefficients for d^2, d^1, and d^0 are -2, 2, and 4.
q_coeffs = [-2, 2, 4]

# The problem asks for the sum of the absolute values of these integer coefficients.

# First, let's print the equation as requested.
print("The difference is expressed as the fraction p(d)/q(d):")
print(f"p(d) = ({p_coeffs[0]})d^2 + ({p_coeffs[1]})d + ({p_coeffs[2]})")
print(f"q(d) = ({q_coeffs[0]})d^2 + ({q_coeffs[1]})d + ({q_coeffs[2]})")
print("-" * 20)

# Calculate the sum of the absolute values of the coefficients for p(d).
sum_abs_p = sum(abs(c) for c in p_coeffs)

# Calculate the sum of the absolute values of the coefficients for q(d).
sum_abs_q = sum(abs(c) for c in q_coeffs)

# The total sum is the sum of these two results.
total_sum = sum_abs_p + sum_abs_q

print(f"The coefficients of the numerator p(d) are: {p_coeffs}")
print(f"The sum of their absolute values is: {sum_abs_p}")

print(f"The coefficients of the denominator q(d) are: {q_coeffs}")
print(f"The sum of their absolute values is: {sum_abs_q}")
print("-" * 20)

print(f"The total sum of the absolute values of all coefficients is: {total_sum}")
