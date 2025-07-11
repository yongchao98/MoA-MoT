import sympy as sp

# Define t and theta as a function of t
t = sp.Symbol('t')
theta = sp.Function('theta')(t)

# Define f as a function of t
f = sp.Function('f')(t)

# Calculate the derivative of f
f_prime = f.diff(t)

# Expression for theta'(t) as found in the derivation
theta_prime_expr = f * (sp.cos(theta)**2) + (f_prime / f) * sp.cos(theta) * sp.sin(theta)

# The question asks for the value of theta'(t). 
# We present the derived formula as the answer.
# Let's print the terms to match the answer choices.
term1_str = f"f(t)*cos^2(theta(t))"
term2_str = f"(f'(t)/f(t))*cos(theta(t))*sin(theta(t))"
final_equation = f"{term1_str} + {term2_str}"

print("The derivative of theta with respect to t, theta'(t), is given by the equation:")
print(f"theta'(t) = {final_equation}")

# We can also pretty print the symbolic expression from sympy
print("\nSymbolic expression from sympy:")
sp.pprint(theta_prime_expr)