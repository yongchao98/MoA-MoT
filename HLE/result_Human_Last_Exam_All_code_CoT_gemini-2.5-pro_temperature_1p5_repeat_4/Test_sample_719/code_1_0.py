import sympy as sp

# Define symbols
t = sp.Symbol('t')
f = sp.Function('f')(t)
theta = sp.Function('theta')(t)
r = sp.Function('r')(t)

# Define components a(t) and b(t)
a = r * sp.cos(theta)
b = r * sp.sin(theta)

# Differential equations for a and b from the Lie derivative calculation
f_prime = f.diff(t)
# dadt = -a * f'/f
# dbdt = a * f
dadt_expr = -a * f_prime / f
dbdt_expr = a * f

# Formula for the derivative of the angle theta
theta_prime_formula = (a * dbdt_expr - b * dadt_expr) / (a**2 + b**2)

print("The formula for the time derivative of the angle theta is:")
print("theta'(t) = (a(t) * b'(t) - b(t) * a'(t)) / (a(t)**2 + b(t)**2)\n")

print("From the geometry of the linearized flow, we derive the following ODEs for a(t) and b(t):")
print(f"a'(t) = {sp.printing.latex(dadt_expr.subs(r * sp.cos(theta), sp.Symbol('a')))}")
print(f"b'(t) = {sp.printing.latex(dbdt_expr.subs(r * sp.cos(theta), sp.Symbol('a')))}\n")

print("Substituting the ODEs into the formula for theta'(t):")
# Substitute dadt and dbdt expressions
theta_prime_intermediate = (a * (a * f) - b * (-a * f_prime / f)) / (a**2 + b**2)
print("theta'(t) = (a(t) * (a(t)*f(t)) - b(t) * (-a(t)*f'(t)/f(t))) / (a(t)**2 + b(t)**2)")
print("theta'(t) = (a(t)**2 * f(t) + a(t)*b(t)*f'(t)/f(t)) / (a(t)**2 + b(t)**2)\n")

print("Now, substitute a(t) = r(t)*cos(theta(t)) and b(t) = r(t)*sin(theta(t)):")
# Substitute a and b with their polar forms
final_expr_num = (r**2 * sp.cos(theta)**2 * f + r**2 * sp.cos(theta) * sp.sin(theta) * f_prime / f)
final_expr_den = r**2
print("theta'(t) = (r(t)**2*cos(theta(t))**2*f(t) + r(t)**2*cos(theta(t))*sin(theta(t))*f'(t)/f(t)) / r(t)**2\n")

# Simplify by canceling r**2
final_expr = final_expr_num / final_expr_den
print("After simplification, we get the final result:")
print(f"theta'(t) = {sp.printing.pretty(final_expr)}")

# To be explicit about the final equation, print each term
term1 = f * sp.cos(theta)**2
term2 = (f_prime / f) * sp.cos(theta) * sp.sin(theta)
print("\nThe final equation is:")
print(f"theta'(t) = {sp.printing.pretty(term1)} + {sp.printing.pretty(term2)}")