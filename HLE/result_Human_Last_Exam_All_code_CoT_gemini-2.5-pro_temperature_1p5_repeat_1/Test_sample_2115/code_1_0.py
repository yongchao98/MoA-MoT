import sympy
from sympy import exp, log, E

# Based on the analysis, the integral's value is conserved over time.
# We calculate the integral at t=0.
# The integral is: I = -3 * Integral from 0 to 1 of (2*e^(2x) + e^x) / (e^(2x) + e^x + 1) dx
# After substitution w = exp(x), this becomes:
# I = -3 * Integral from 1 to e of (2w + 1) / (w^2 + w + 1) dw
# The antiderivative is -3 * ln(w^2 + w + 1).

# Define the symbolic variables
w = sympy.Symbol('w')
x = sympy.Symbol('x')

# Define the antiderivative in terms of w
antiderivative_w = -3 * log(w**2 + w + 1)

# Evaluate the definite integral using the limits w=1 and w=E
val_upper = antiderivative_w.subs(w, E)
val_lower = antiderivative_w.subs(w, 1)
integral_result = val_upper - val_lower

# Numerically evaluate the result
numerical_result = integral_result.evalf()

# Print the components of the final equation
c = -3
term1_expr = E**2 + E + 1
term2_val = 3

print("The problem reduces to calculating the integral at t=0, which is of the form C * (ln(A) - ln(B)).")
print(f"The equation for the result is: {c} * (ln({term1_expr}) - ln({term2_val}))")
print("\nWhere:")
print(f"C = {c}")
print(f"A = e^2 + e + 1 = {term1_expr.evalf()}")
print(f"B = {term2_val}")
print("\nThis simplifies to the exact form:")
print(f"{c} * ln((e^2 + e + 1)/3)")
print("\nFinal numerical result:")
print(numerical_result)