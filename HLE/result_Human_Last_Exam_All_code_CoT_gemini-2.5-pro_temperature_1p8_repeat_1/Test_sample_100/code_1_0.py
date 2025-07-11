import sympy

# The problem simplifies to calculating the sum of four definite integrals.
# We will calculate each of these values and sum them up.

# Integral 1: integral from 0 to oo of p/(exp(p)-1) dp
# This is Gamma(2)*zeta(2)
pi = sympy.pi
val1 = sympy.gamma(2) * sympy.zeta(2)

# Integral 2: integral from 0 to oo of p**7/(exp(p)-1) dp
# This is Gamma(8)*zeta(8)
val2 = sympy.gamma(8) * sympy.zeta(8)

# Integral 3: integral from 0 to oo of (exp(p/4)-exp(-p/4))/(2*(exp(p)-1)) dp
# As derived in the explanation, this evaluates to 2 - pi/2
val3 = 2 - pi/2

# Integral 4: integral from 0 to oo of p*exp(-p)/(exp(p)-1) dp
# As derived, this evaluates to Gamma(2)*zeta(2) - Gamma(2)
val4 = val1 - sympy.gamma(2)

# The total integral is the sum of these four values.
total_value = val1 + val2 + val3 + val4

# Simplify the final expression
final_expression = sympy.simplify(total_value)

# The problem requires outputting each number in the final equation.
# The final expression is of the form A*pi^8 + B*pi^2 + C*pi + D
# We will extract the coefficients and print them.

# Expand the expression to make it easier to extract coefficients
expanded_expr = sympy.expand(final_expression)

# Extract coefficients
coeff_pi8 = expanded_expr.coeff(pi**8)
coeff_pi2 = expanded_expr.coeff(pi**2)
coeff_pi = expanded_expr.coeff(pi)
# The constant term D is what remains
const_term = expanded_expr - coeff_pi8*pi**8 - coeff_pi2*pi**2 - coeff_pi*pi

# Print the final equation with coefficients
print("The final result is an expression of the form: A * pi**8 + B * pi**2 + C * pi + D")
print(f"The equation is: ({coeff_pi8}) * pi^8 + ({coeff_pi2}) * pi^2 + ({coeff_pi}) * pi + ({const_term})")
print("This corresponds to the expression:")
print(final_expression)
