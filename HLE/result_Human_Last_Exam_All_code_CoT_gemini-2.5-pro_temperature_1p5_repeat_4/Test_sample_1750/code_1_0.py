import sympy

# Define the symbolic variable for integration
x = sympy.Symbol('x')

# Define the simplified integrand
# The expression is (2*sin(x)*cos(3x))**50
integrand = (2 * sympy.sin(x) * sympy.cos(3*x))**50

# Perform the symbolic integration from 0 to pi
result = sympy.integrate(integrand, (x, 0, sympy.pi))

# The result is a rational number multiplied by pi.
# We extract the numerator and denominator of this rational number.
# result = coefficient * pi
# We can get the coefficient by dividing by pi
coefficient = result / sympy.pi

# sympy represents rational numbers with p and q attributes for numerator and denominator
numerator = coefficient.p
denominator = coefficient.q

# Print the numbers in the final equation as requested
print("The value of the integral is (N * pi) / D where:")
print(f"N = {numerator}")
print(f"D = {denominator}")

# To output the final answer in the required format, let's represent it clearly
final_answer_string = f"{numerator}*pi/{denominator}"
print(f"The integral evaluates to: {final_answer_string}")

# Final Answer for the platform
# The platform expects the final answer string in a specific format.
# Let's construct it. Assuming the fraction containing pi is the desired output.
final_answer = f"<<<{numerator}*pi/{denominator}>>>"
# print(final_answer) # This would be for the final submission.